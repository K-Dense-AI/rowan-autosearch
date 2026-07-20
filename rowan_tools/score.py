"""Submit candidate SMILES to Rowan, score them against the objective, and
append results to the run's iteration record.

The agent calls this once per iteration with a list of candidate SMILES,
their parents, and design notes. We submit in parallel, parse the metric
out of object_data, check constraints (locally via RDKit), and write a new
iterations/NNNN.json file.
"""
from __future__ import annotations

import json
import os
import sys
import threading
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from inspect import signature
from pathlib import Path
from typing import Any, Callable

import click
from dotenv import load_dotenv

from . import state

load_dotenv(Path(__file__).resolve().parent.parent / ".env")


# ---------------------------------------------------------------------------
# Metric extraction
# ---------------------------------------------------------------------------

# For each workflow_type, a list of (metric_name, dotted_path) candidates that
# score.py will try when the user did not pin metric_path in config. Order
# matters — first hit wins. The agent can always override via --metric-path
# in `rowan-state init`.
DEFAULT_METRIC_PATHS: dict[str, list[tuple[str, str]]] = {
    "solubility": [
        ("logS", "solubilities.O.solubilities.0"),
        ("logS", "solubility.logS"),
        ("logS_25C", "solubility.logS_25C"),
        ("logS", "logS"),
    ],
    "pka": [
        ("pka", "strongest_acid"),
        ("pka", "pka_values.0"),
        ("pka", "pka"),
    ],
    "redox_potential": [
        ("redox_potential", "reduction_potential"),
        ("redox_potential", "oxidation_potential"),
        ("redox_potential", "redox_potential"),
    ],
    "admet": [
        ("logP", "properties.logP"),
        ("logP", "logP"),
        ("logD", "properties.logD"),
        ("logD", "logD"),
        ("solubility", "properties.solubility"),
        ("solubility", "solubility"),
    ],
    "bde": [
        ("bde", "bde"),
        ("bde", "weakest_bond.bde"),
    ],
    "descriptors": [
        ("logP", "descriptors.logP"),
        ("tpsa", "descriptors.tpsa"),
    ],
    "docking": [
        ("docking_score", "scores.0.score"),
        ("docking_score", "scores.0"),
    ],
    "batch_docking": [
        ("docking_score", "scores"),
    ],
}


def walk(obj: Any, dotted: str) -> Any:
    """Walk obj by dot-path; supports list indices via integer segments."""
    cur = obj
    for seg in dotted.split("."):
        if cur is None:
            return None
        if isinstance(cur, list):
            try:
                cur = cur[int(seg)]
                continue
            except (ValueError, IndexError):
                return None
        if isinstance(cur, dict):
            cur = cur.get(seg)
            continue
        return None
    return cur


def extract_metric(workflow_type: str, primary_metric: str,
                   metric_path: str | None,
                   object_data: dict[str, Any]) -> tuple[float | None, str | None]:
    """Return (value, path_used). Tries pinned path first, then heuristics."""
    if metric_path:
        v = walk(object_data, metric_path)
        if isinstance(v, (int, float)):
            return float(v), metric_path

    for metric_name, path in DEFAULT_METRIC_PATHS.get(workflow_type, []):
        if metric_name != primary_metric:
            continue
        v = walk(object_data, path)
        if isinstance(v, (int, float)):
            return float(v), path

    # Last-resort: top-level key match
    v = object_data.get(primary_metric)
    if isinstance(v, (int, float)):
        return float(v), primary_metric

    return None, None


def _term_value(term: dict[str, Any], metrics: dict[str, float],
                object_data: dict[str, Any]) -> tuple[float | None, str | None]:
    """Resolve an optimization term from local metrics or Rowan object_data."""
    source = term.get("source", "rowan")
    metric = term.get("metric") or term.get("name")
    path = term.get("path")

    if source == "local":
        value = metrics.get(metric) if metric else None
        return (float(value), metric) if isinstance(value, (int, float)) else (None, metric)

    if path:
        value = walk(object_data, path)
        return (float(value), path) if isinstance(value, (int, float)) else (None, path)

    if metric:
        value = object_data.get(metric)
        if isinstance(value, (int, float)):
            return float(value), metric
        value = metrics.get(metric)
        if isinstance(value, (int, float)):
            return float(value), metric

    return None, path or metric


def evaluate_optimization_objective(
    cfg: dict[str, Any],
    metrics: dict[str, float],
    object_data: dict[str, Any],
) -> tuple[float | None, dict[str, Any] | None, str | None]:
    """Evaluate optional weighted objective terms.

    Config shape:
      "optimization_objective": {
        "name": "objective",
        "terms": [
          {"name": "logS", "path": "solubilities.O.solubilities.0", "goal": "maximize", "weight": 1.0},
          {"name": "mw", "source": "local", "metric": "mw", "goal": "minimize", "weight": 0.01},
          {"name": "pKa_target", "path": "strongest_acid", "target": 4.5, "weight": 1.0}
        ]
      }
    """
    objective = cfg.get("optimization_objective") or {}
    terms = objective.get("terms") or cfg.get("objective_terms") or []
    if not terms:
        return None, None, None

    total = 0.0
    details: dict[str, Any] = {}
    for i, term in enumerate(terms):
        if not isinstance(term, dict):
            return None, None, f"objective term {i} is not an object"

        value, path_used = _term_value(term, metrics, object_data)
        name = str(term.get("name") or term.get("metric") or path_used or f"term_{i + 1}")
        if value is None:
            return None, details, f"objective term {name!r} could not be resolved from {path_used!r}"

        weight = float(term.get("weight", 1.0))
        goal = term.get("goal") or term.get("direction") or "maximize"
        if term.get("target") is not None or goal == "target":
            if term.get("target") is None:
                return None, details, f"objective term {name!r} has goal 'target' but no target value"
            target = float(term["target"])
            contribution = -weight * abs(value - target)
        elif goal == "minimize":
            contribution = -weight * value
        else:
            contribution = weight * value

        total += contribution
        details[name] = {
            "value": value,
            "path": path_used,
            "goal": goal,
            "weight": weight,
            "contribution": contribution,
        }

    return total, details, None


# ---------------------------------------------------------------------------
# Local descriptors (constraint checking, doesn't burn Rowan credits)
# ---------------------------------------------------------------------------

def local_descriptors(smiles: str) -> dict[str, float]:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp_crippen": float(Crippen.MolLogP(mol)),
        "tpsa": float(Descriptors.TPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotb": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "rings": int(mol.GetRingInfo().NumRings()),
    }


CONSTRAINT_RULES: dict[str, tuple[str, str]] = {
    # constraint_key: (descriptor_key, comparator)  comparator in {"<=", ">="}
    "mw_max": ("mw", "<="),
    "mw_min": ("mw", ">="),
    "logp_max": ("logp_crippen", "<="),
    "logp_min": ("logp_crippen", ">="),
    "tpsa_max": ("tpsa", "<="),
    "tpsa_min": ("tpsa", ">="),
    "hbd_max": ("hbd", "<="),
    "hba_max": ("hba", "<="),
    "rotb_max": ("rotb", "<="),
    "ha_max": ("heavy_atoms", "<="),
    "ha_min": ("heavy_atoms", ">="),
    "rings_max": ("rings", "<="),
}


def check_constraints(desc: dict[str, float],
                      constraints: dict[str, Any]) -> tuple[bool, list[str]]:
    failed = []
    for k, target in constraints.items():
        rule = CONSTRAINT_RULES.get(k)
        if not rule:
            continue
        dkey, op = rule
        v = desc.get(dkey)
        if v is None:
            continue
        if op == "<=" and v > float(target):
            failed.append(f"{k}: {dkey}={v:.2f} > {target}")
        if op == ">=" and v < float(target):
            failed.append(f"{k}: {dkey}={v:.2f} < {target}")
    return (len(failed) == 0), failed


def canonical_smiles(smiles: str) -> str | None:
    """Canonical SMILES used for hard dedup at the scoring boundary."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def submitted_history(run_id: str) -> dict[str, dict[str, Any]]:
    """Canonical SMILES -> latest record that has already consumed Rowan work."""
    seen: dict[str, dict[str, Any]] = {}
    for candidate in state.all_candidates(run_id):
        if candidate.get("score") is None and not candidate.get("workflow_uuid"):
            continue
        key = canonical_smiles(candidate.get("smiles", ""))
        if key:
            seen[key] = candidate
    return seen


def preflight_candidates(
    run_id: str,
    candidates: list[dict[str, str]],
    cfg: dict[str, Any],
    *,
    force_rescore: bool = False,
    score_constraint_failures: bool = False,
) -> tuple[list[dict[str, Any]], list[int]]:
    """Apply hard local gates and return candidate stubs plus submit indexes."""
    history = submitted_history(run_id)
    batch_seen: dict[str, int] = {}
    records: list[dict[str, Any]] = []
    submit_indexes: list[int] = []

    for candidate in candidates:
        smiles = candidate["smiles"]
        parent = candidate.get("parent_smiles") or None
        note = candidate.get("design_note", "")
        desc = local_descriptors(smiles)
        record: dict[str, Any] = {
            "smiles": smiles,
            "parent_smiles": parent,
            "design_note": note,
            "prepared_at": state.now_iso(),
            "metrics": dict(desc),
            "score": None,
            "status": "pending",
        }

        key = canonical_smiles(smiles)
        if not desc or key is None:
            record.update(
                status="skipped",
                error="RDKit failed to parse SMILES",
                satisfies_constraints=False,
                constraint_failures=["invalid SMILES"],
                skip_reason="invalid SMILES",
            )
            records.append(record)
            continue

        record["canonical_smiles"] = key
        ok, failed = check_constraints(desc, cfg.get("constraints", {}))
        record["satisfies_constraints"] = ok
        record["constraint_failures"] = failed

        if not ok and not score_constraint_failures:
            record.update(
                status="skipped",
                skip_reason="failed local constraints",
            )
            records.append(record)
            continue

        if key in batch_seen:
            record.update(
                status="skipped",
                skip_reason="duplicate candidate in this batch",
                duplicate_of_index=batch_seen[key],
            )
            records.append(record)
            continue

        prior = history.get(key)
        if prior is not None and not force_rescore:
            record.update(
                status="skipped",
                skip_reason="already submitted to Rowan",
                prior_iter=prior.get("iter"),
                prior_score=prior.get("score"),
                prior_workflow_uuid=prior.get("workflow_uuid"),
            )
            records.append(record)
            continue

        batch_seen[key] = len(records)
        submit_indexes.append(len(records))
        records.append(record)

    return records, submit_indexes


# ---------------------------------------------------------------------------
# Rowan submission
# ---------------------------------------------------------------------------

def _workflow_params(cfg: dict[str, Any]) -> dict[str, Any]:
    return dict(cfg.get("workflow_params") or cfg.get("rowan_workflow_kwargs") or {})


def _submit_with_smiles(submit: Any, smiles: str, name: str,
                        params: dict[str, Any]) -> Any:
    """Call a Rowan submitter using the molecule parameter name it exposes."""
    sig = signature(submit)
    kwargs = dict(params)
    if "name" in sig.parameters:
        kwargs.setdefault("name", name)

    for molecule_param in ("initial_smiles", "initial_molecule", "molecule"):
        if molecule_param in sig.parameters:
            kwargs.setdefault(molecule_param, smiles)
            return submit(**kwargs)

    # Fallback for future helpers that still accept the molecule positionally.
    return submit(smiles, **kwargs)


class RowanWorkflowError(RuntimeError):
    """Failure after Rowan may have accepted a workflow."""

    def __init__(self, message: str, workflow_uuid: str | None = None) -> None:
        super().__init__(message)
        self.workflow_uuid = workflow_uuid


def submit_one(
    smiles: str,
    cfg: dict[str, Any],
    name: str,
    *,
    workflow_uuid: str | None = None,
    on_submitted: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Submit or retrieve a Rowan workflow, checkpoint its UUID, and wait."""
    import rowan

    workflow_type = cfg["workflow_type"]
    params = _workflow_params(cfg)
    if workflow_uuid:
        workflow = rowan.retrieve_workflow(workflow_uuid)
    else:
        if workflow_type == "solubility":
            defaults = {"method": "kingfisher", "solvents": ["water"], "temperatures": [298.15]}
            defaults.update(params)
            workflow = _submit_with_smiles(rowan.submit_solubility_workflow, smiles, name, defaults)
        else:
            submit = getattr(rowan, f"submit_{workflow_type}_workflow", None)
            if submit is None:
                raise AttributeError(f"rowan has no submit function for workflow {workflow_type!r}")
            workflow = _submit_with_smiles(submit, smiles, name, params)

    uuid = str(workflow.uuid)
    try:
        if on_submitted is not None:
            on_submitted(uuid)
        result = workflow.result(wait=True, poll_interval=5)
    except Exception as e:
        raise RowanWorkflowError(str(e), uuid) from e
    return {
        "uuid": uuid,
        "object_data": result.data,
        "workflow_type": workflow_type,
    }


def score_candidate(cand: dict[str, str], cfg: dict[str, Any],
                    iter_n: int, *,
                    workflow_uuid: str | None = None,
                    on_submitted: Callable[[str], None] | None = None,
                    score_constraint_failures: bool = False) -> dict[str, Any]:
    smiles = cand["smiles"]
    parent = cand.get("parent_smiles") or None
    note = cand.get("design_note", "")

    desc = local_descriptors(smiles)
    if not desc:
        return {
            "smiles": smiles, "parent_smiles": parent, "design_note": note,
            "error": "RDKit failed to parse SMILES", "score": None,
            "metrics": {}, "satisfies_constraints": False,
            "constraint_failures": ["invalid SMILES"],
            "status": "skipped", "skip_reason": "invalid SMILES",
        }

    ok, failed = check_constraints(desc, cfg.get("constraints", {}))

    out: dict[str, Any] = {
        "smiles": smiles,
        "parent_smiles": parent,
        "design_note": note,
        "scored_at": state.now_iso(),
        "metrics": dict(desc),
        "satisfies_constraints": ok,
        "constraint_failures": failed,
        "score": None,
        "status": "pending",
    }

    if not ok and not score_constraint_failures:
        out["status"] = "skipped"
        out["skip_reason"] = "failed local constraints"
        return out

    if workflow_uuid:
        out["workflow_uuid"] = workflow_uuid
    name = f"{cfg['run_id']}_iter{iter_n:04d}_{smiles[:24]}"
    try:
        if workflow_uuid is None and on_submitted is None:
            result = submit_one(smiles, cfg, name)
        else:
            result = submit_one(
                smiles,
                cfg,
                name,
                workflow_uuid=workflow_uuid,
                on_submitted=on_submitted,
            )
    except Exception as e:
        out["error"] = f"Rowan submit failed: {e}"
        out["traceback"] = traceback.format_exc()
        recovered_uuid = getattr(e, "workflow_uuid", None)
        if recovered_uuid:
            out["workflow_uuid"] = recovered_uuid
        out["status"] = "failed"
        return out

    obj = result.get("object_data", result) if isinstance(result, dict) else {}
    out["workflow_uuid"] = result.get("uuid") if isinstance(result, dict) else None
    out["object_data"] = obj  # keep full payload for agent inspection

    val, path = extract_metric(
        cfg["workflow_type"],
        cfg["primary_metric"],
        cfg.get("metric_path"),
        obj,
    )
    out["metric_path_used"] = path
    if val is not None:
        out["metrics"][cfg["primary_metric"]] = val
        out["score"] = val

    objective_score, objective_details, objective_error = evaluate_optimization_objective(
        cfg, out["metrics"], obj
    )
    if objective_details is not None:
        out["optimization_objective"] = objective_details
    if objective_score is not None:
        objective_name = (cfg.get("optimization_objective") or {}).get("name", "objective_score")
        out["metrics"][objective_name] = objective_score
        out["score"] = objective_score
    elif objective_error:
        out["error"] = f"Could not evaluate optimization objective: {objective_error}"
    elif val is None:
        out["error"] = (
            f"Could not extract metric {cfg['primary_metric']!r} from object_data. "
            f"Inspect candidates[*].object_data and set metric_path in config.json."
        )

    out["status"] = "complete" if out.get("score") is not None else "failed"
    return out


def score_batch_docking(
    candidates: list[dict[str, str]],
    cfg: dict[str, Any],
    iter_n: int,
    *,
    workflow_uuid: str | None = None,
    on_submitted: Callable[[str], None] | None = None,
    score_constraint_failures: bool = False,
) -> list[dict[str, Any]]:
    """Score a ligand set with Rowan's batch_docking workflow."""
    import rowan

    scored = []
    smiles_list = []
    submitted_indexes = []
    for i, cand in enumerate(candidates):
        smiles = cand["smiles"]
        desc = local_descriptors(smiles)
        parent = cand.get("parent_smiles") or None
        note = cand.get("design_note", "")
        if not desc:
            scored.append({
                "smiles": smiles, "parent_smiles": parent, "design_note": note,
                "error": "RDKit failed to parse SMILES", "score": None,
                "metrics": {}, "satisfies_constraints": False,
                "constraint_failures": ["invalid SMILES"],
                "status": "skipped", "skip_reason": "invalid SMILES",
            })
            continue

        ok, failed = check_constraints(desc, cfg.get("constraints", {}))
        out = {
            "smiles": smiles,
            "parent_smiles": parent,
            "design_note": note,
            "scored_at": state.now_iso(),
            "metrics": dict(desc),
            "satisfies_constraints": ok,
            "constraint_failures": failed,
            "score": None,
            "status": "pending",
        }
        scored.append(out)
        if not ok and not score_constraint_failures:
            out["status"] = "skipped"
            out["skip_reason"] = "failed local constraints"
            continue
        smiles_list.append(smiles)
        submitted_indexes.append(i)

    if not smiles_list:
        return scored

    params = _workflow_params(cfg)
    name = f"{cfg['run_id']}_iter{iter_n:04d}_batch_docking"
    workflow = None
    try:
        if workflow_uuid:
            workflow = rowan.retrieve_workflow(workflow_uuid)
        else:
            workflow = rowan.submit_batch_docking_workflow(
                smiles_list=smiles_list,
                name=name,
                **params,
            )
        uuid = str(workflow.uuid)
        for i in submitted_indexes:
            scored[i]["workflow_uuid"] = uuid
            scored[i]["status"] = "submitted"
        if on_submitted is not None:
            on_submitted(uuid)
        result = workflow.result(wait=True, poll_interval=5)
        object_data = result.data
        batch_scores = result.scores
    except Exception as e:
        error = f"Rowan batch_docking submit failed: {e}"
        tb = traceback.format_exc()
        for i in submitted_indexes:
            scored[i]["error"] = error
            scored[i]["traceback"] = tb
            scored[i]["status"] = "failed"
            if workflow is not None:
                scored[i]["workflow_uuid"] = str(workflow.uuid)
        return scored

    canonical_scores = {
        canonical_smiles(returned_smiles): value
        for returned_smiles, value in batch_scores.items()
        if canonical_smiles(returned_smiles) is not None
    }
    for i in submitted_indexes:
        out = scored[i]
        smiles = out["smiles"]
        val = batch_scores.get(smiles)
        if val is None:
            val = canonical_scores.get(canonical_smiles(smiles))
        out["workflow_uuid"] = str(workflow.uuid)
        out["object_data"] = object_data
        out["metric_path_used"] = f"scores.{smiles}"
        if isinstance(val, (int, float)):
            out["metrics"][cfg["primary_metric"]] = float(val)
            out["score"] = float(val)
        else:
            out["error"] = f"No batch docking score returned for {smiles}"

        objective_score, objective_details, objective_error = evaluate_optimization_objective(
            cfg, out["metrics"], object_data
        )
        if objective_details is not None:
            out["optimization_objective"] = objective_details
        if objective_score is not None:
            objective_name = (cfg.get("optimization_objective") or {}).get("name", "objective_score")
            out["metrics"][objective_name] = objective_score
            out["score"] = objective_score
        elif objective_error:
            out["error"] = f"Could not evaluate optimization objective: {objective_error}"

        out["status"] = "complete" if out.get("score") is not None else "failed"

    return scored


def rebuild_report(run_id: str) -> Path:
    """Regenerate the HTML report after an iteration is written."""
    from viz.build_report import build_report

    return build_report(run_id)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _iteration_best(
    candidates: list[dict[str, Any]],
    direction: str,
) -> dict[str, Any] | None:
    valid = [
        candidate for candidate in candidates
        if candidate.get("score") is not None
        and candidate.get("satisfies_constraints", True)
    ]
    if not valid:
        return None
    chooser = max if direction == "maximize" else min
    return chooser(valid, key=lambda candidate: candidate["score"])


def _incomplete_iterations(run_id: str) -> list[int]:
    return [
        n for n in state.list_iterations(run_id)
        if state.load_iter(run_id, n).get("status") in {"in_progress", "incomplete"}
    ]


def _print_candidate_result(
    candidate: dict[str, Any],
    index: int,
    total: int,
    metric: str,
) -> None:
    if candidate.get("score") is not None:
        result = f"{metric}={candidate['score']:.3f}"
    elif candidate.get("status") == "skipped":
        result = f"SKIPPED: {candidate.get('skip_reason', 'local gate')}"
    else:
        result = "FAILED"
    click.echo(f"  [{index + 1}/{total}] {candidate['smiles'][:40]:40s}  {result}")


def _finalize_iteration(
    run_id: str,
    cfg: dict[str, Any],
    record: dict[str, Any],
) -> None:
    best = _iteration_best(record["candidates"], cfg["objective_direction"])
    record["iteration_best"] = (
        {"smiles": best["smiles"], "score": best["score"]} if best else None
    )
    unfinished = any(
        candidate.get("status") not in {"complete", "skipped"}
        for candidate in record["candidates"]
    )
    record["status"] = "incomplete" if unfinished else "complete"
    record["updated_at"] = state.now_iso()
    if not unfinished:
        record["completed_at"] = record["updated_at"]
    state.save_iter(run_id, record["iter"], record)


def _score_iteration(
    run_id: str,
    cfg: dict[str, Any],
    record: dict[str, Any],
    submit_indexes: list[int],
    *,
    max_workers: int,
    score_constraint_failures: bool,
) -> None:
    """Score eligible candidates and checkpoint UUIDs/results as they arrive."""
    iter_n = record["iter"]
    checkpoint_lock = threading.Lock()

    def checkpoint_submitted(indexes: list[int], uuid: str) -> None:
        with checkpoint_lock:
            for index in indexes:
                candidate = record["candidates"][index]
                candidate["workflow_uuid"] = uuid
                candidate["status"] = "submitted"
                candidate["submitted_at"] = state.now_iso()
            record["updated_at"] = state.now_iso()
            state.save_iter(run_id, iter_n, record)

    if cfg["workflow_type"] == "batch_docking" and submit_indexes:
        active = [record["candidates"][index] for index in submit_indexes]
        uuids = {
            candidate["workflow_uuid"]
            for candidate in active
            if candidate.get("workflow_uuid")
        }
        if len(uuids) > 1:
            raise click.ClickException(
                "Cannot resume batch docking: candidates reference multiple workflow UUIDs."
            )
        workflow_uuid = next(iter(uuids), None)
        scored = score_batch_docking(
            active,
            cfg,
            iter_n,
            workflow_uuid=workflow_uuid,
            on_submitted=lambda uuid: checkpoint_submitted(submit_indexes, uuid),
            score_constraint_failures=score_constraint_failures,
        )
        with checkpoint_lock:
            for index, result in zip(submit_indexes, scored):
                merged = dict(record["candidates"][index])
                merged.update(result)
                record["candidates"][index] = merged
            record["updated_at"] = state.now_iso()
            state.save_iter(run_id, iter_n, record)
        for index in submit_indexes:
            _print_candidate_result(
                record["candidates"][index],
                index,
                len(record["candidates"]),
                cfg["primary_metric"],
            )
    elif submit_indexes:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for index in submit_indexes:
                candidate = record["candidates"][index]
                future = executor.submit(
                    score_candidate,
                    candidate,
                    cfg,
                    iter_n,
                    workflow_uuid=candidate.get("workflow_uuid"),
                    on_submitted=lambda uuid, i=index: checkpoint_submitted([i], uuid),
                    score_constraint_failures=score_constraint_failures,
                )
                futures[future] = index

            for future in as_completed(futures):
                index = futures[future]
                try:
                    result = future.result()
                except Exception as e:
                    result = {
                        **record["candidates"][index],
                        "error": str(e),
                        "traceback": traceback.format_exc(),
                        "score": None,
                        "status": "failed",
                    }
                with checkpoint_lock:
                    merged = dict(record["candidates"][index])
                    merged.update(result)
                    record["candidates"][index] = merged
                    record["updated_at"] = state.now_iso()
                    state.save_iter(run_id, iter_n, record)
                _print_candidate_result(
                    record["candidates"][index],
                    index,
                    len(record["candidates"]),
                    cfg["primary_metric"],
                )

    _finalize_iteration(run_id, cfg, record)


def _rebuild_report_with_warning(run_id: str) -> None:
    try:
        report_path = rebuild_report(run_id)
    except Exception as e:
        click.echo(f"\nWarning: could not regenerate report.html: {e}", err=True)
        click.echo(
            f"Run `rowan-report --run {run_id}` after fixing the report error.",
            err=True,
        )
    else:
        click.echo(f"\nUpdated report: {report_path}")
        click.echo(f"Open it: file://{report_path.resolve()}")


@click.command()
@click.option("--run", "run_id", required=True)
@click.option("--rationale", required=True,
              help="Why these candidates? (Saved with the iteration.)")
@click.option("--candidate", "candidates", multiple=True,
              help="Candidate as 'SMILES|PARENT_SMILES|DESIGN_NOTE'. "
                   "PARENT_SMILES may be empty. Repeat once per candidate.")
@click.option("--decision", default="continue",
              type=click.Choice(["continue", "converged", "stuck"]),
              help="Agent's decision after seeing scores. Default 'continue' — "
                   "update later by re-running with --decision.")
@click.option("--max-workers", default=4, type=int,
              help="Parallel Rowan submissions.")
@click.option("--dry-run", is_flag=True, help="Validate inputs without calling Rowan.")
@click.option("--resume", is_flag=True,
              help="Resume the latest incomplete iteration using saved Rowan UUIDs.")
@click.option("--force-rescore", is_flag=True,
              help="Allow a molecule already submitted in this run to be scored again.")
@click.option("--score-constraint-failures", is_flag=True,
              help="Explicitly submit candidates that fail local constraints.")
@click.option("--allow-over-budget", is_flag=True,
              help="Allow a new scoring iteration beyond config max_iterations.")
def cli(run_id, rationale, candidates, decision, max_workers, dry_run, resume,
        force_rescore, score_constraint_failures, allow_over_budget):
    """Score candidates and append a new iteration to the run."""
    cfg = state.load_config(run_id)

    parsed = []
    for c in candidates:
        parts = c.split("|", 2)
        smiles = parts[0].strip()
        parent = parts[1].strip() if len(parts) > 1 else ""
        note = parts[2].strip() if len(parts) > 2 else ""
        if not smiles:
            raise click.ClickException(f"Empty SMILES in --candidate {c!r}")
        parsed.append({"smiles": smiles, "parent_smiles": parent, "design_note": note})

    if resume:
        if parsed:
            raise click.ClickException("--resume cannot be combined with --candidate.")
        if dry_run:
            raise click.ClickException("--resume cannot be combined with --dry-run.")
        incomplete = _incomplete_iterations(run_id)
        if not incomplete:
            raise click.ClickException("No incomplete iteration is available to resume.")
        iter_n = incomplete[-1]
        iter_record = state.load_iter(run_id, iter_n)
        submit_indexes = [
            index for index, candidate in enumerate(iter_record.get("candidates", []))
            if candidate.get("status") not in {"complete", "skipped"}
        ]
        if not submit_indexes:
            raise click.ClickException(
                f"Iteration {iter_n} has no pending or failed candidates to resume."
            )
        if not os.environ.get("ROWAN_API_KEY"):
            raise click.ClickException(
                "ROWAN_API_KEY not set. Copy .env.example to .env and fill it in."
            )
        iter_record.setdefault("resume_events", []).append({
            "at": state.now_iso(),
            "rationale": rationale,
        })
        iter_record["decision"] = decision
        if score_constraint_failures:
            iter_record["score_constraint_failures"] = True
        iter_record["status"] = "in_progress"
        iter_record["updated_at"] = state.now_iso()
        state.save_iter(run_id, iter_n, iter_record)
        click.echo(
            f"[run {run_id}] resuming iter {iter_n}: "
            f"{len(submit_indexes)} candidate(s)"
        )
        _score_iteration(
            run_id,
            cfg,
            iter_record,
            submit_indexes,
            max_workers=max_workers,
            score_constraint_failures=bool(
                iter_record.get("score_constraint_failures", False)
            ),
        )
        click.echo(f"\nSaved {state.iter_path(run_id, iter_n)}")
        best = state.best_so_far(run_id)
        if best:
            click.echo(
                f"Best so far: {best['smiles']}  "
                f"{cfg['primary_metric']}={best['score']:.3f}  (iter {best['iter']})"
            )
        _rebuild_report_with_warning(run_id)
        return

    if not parsed:
        if dry_run:
            click.echo("(dry-run; no candidates to score)")
            return

        existing_iters = state.list_iterations(run_id)
        if decision != "continue":
            if not existing_iters:
                raise click.ClickException("Cannot record a decision before any iterations exist.")
            latest_n = existing_iters[-1]
            latest = state.load_iter(run_id, latest_n)
            latest["decision"] = decision
            latest["decision_rationale"] = rationale
            latest["decision_at"] = state.now_iso()
            state.save_iter(run_id, latest_n, latest)
            click.echo(f"Updated decision for iteration {latest_n} to {decision}.")
            _rebuild_report_with_warning(run_id)
            return

        if not existing_iters:
            raise click.ClickException(
                "Cannot append a wrap-up before any scoring iteration exists."
            )
        wrapup_decision = state.load_iter(
            run_id, existing_iters[-1]
        ).get("decision", "continue")

        iter_record = {
            "created_at": state.now_iso(),
            "rationale": rationale,
            "candidates": [],
            "iteration_best": None,
            "decision": wrapup_decision,
            "status": "complete",
        }
        iter_n = state.reserve_iter(run_id, iter_record)
        click.echo(f"Saved wrap-up iteration {state.iter_path(run_id, iter_n)}")
        _rebuild_report_with_warning(run_id)
        return

    candidate_records, submit_indexes = preflight_candidates(
        run_id,
        parsed,
        cfg,
        force_rescore=force_rescore,
        score_constraint_failures=score_constraint_failures,
    )

    if dry_run:
        click.echo(f"[run {run_id}] preflight: {len(parsed)} candidates")
        for index, candidate in enumerate(candidate_records):
            if index in submit_indexes:
                click.echo(f"  - {candidate['smiles']}  READY")
            else:
                click.echo(
                    f"  - {candidate['smiles']}  "
                    f"SKIPPED: {candidate.get('skip_reason', 'local gate')}"
                )
        click.echo("(dry-run; not calling Rowan)")
        click.echo("(no iteration written)")
        return

    incomplete = _incomplete_iterations(run_id)
    if incomplete:
        raise click.ClickException(
            f"Iteration {incomplete[-1]} is incomplete. Resume it with --resume "
            "before starting a new iteration."
        )
    existing_iters = state.list_iterations(run_id)
    max_iterations = int(cfg.get("max_iterations", 50))
    if len(existing_iters) >= max_iterations and not allow_over_budget:
        raise click.ClickException(
            f"Run has reached max_iterations={max_iterations}. "
            "Use --allow-over-budget to explicitly start another scoring iteration."
        )
    if submit_indexes and not os.environ.get("ROWAN_API_KEY"):
        raise click.ClickException(
            "ROWAN_API_KEY not set. Copy .env.example to .env and fill it in."
        )

    iter_record = {
        "created_at": state.now_iso(),
        "updated_at": state.now_iso(),
        "rationale": rationale,
        "candidates": candidate_records,
        "iteration_best": None,
        "decision": decision,
        "status": "in_progress",
        "force_rescore": force_rescore,
        "score_constraint_failures": score_constraint_failures,
    }
    iter_n = state.reserve_iter(run_id, iter_record)
    iter_record["iter"] = iter_n
    click.echo(
        f"[run {run_id}] iter {iter_n}: {len(parsed)} candidates, "
        f"{len(submit_indexes)} eligible for Rowan"
    )
    for index, candidate in enumerate(candidate_records):
        if index not in submit_indexes:
            _print_candidate_result(
                candidate, index, len(candidate_records), cfg["primary_metric"]
            )
    _score_iteration(
        run_id,
        cfg,
        iter_record,
        submit_indexes,
        max_workers=max_workers,
        score_constraint_failures=score_constraint_failures,
    )
    click.echo(f"\nSaved {state.iter_path(run_id, iter_n)}")

    best = state.best_so_far(run_id)
    if best:
        click.echo(f"Best so far: {best['smiles']}  "
                   f"{cfg['primary_metric']}={best['score']:.3f}  (iter {best['iter']})")

    _rebuild_report_with_warning(run_id)


if __name__ == "__main__":
    cli()
