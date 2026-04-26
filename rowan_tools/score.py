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
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from inspect import signature
from pathlib import Path
from typing import Any

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


def submit_one(smiles: str, cfg: dict[str, Any], name: str) -> dict[str, Any]:
    """Submit a Rowan workflow and wait for its typed result data."""
    import rowan

    workflow_type = cfg["workflow_type"]
    params = _workflow_params(cfg)
    if workflow_type == "solubility":
        defaults = {"method": "kingfisher", "solvents": ["water"], "temperatures": [298.15]}
        defaults.update(params)
        workflow = _submit_with_smiles(rowan.submit_solubility_workflow, smiles, name, defaults)
    else:
        submit = getattr(rowan, f"submit_{workflow_type}_workflow", None)
        if submit is None:
            raise AttributeError(f"rowan has no submit function for workflow {workflow_type!r}")
        workflow = _submit_with_smiles(submit, smiles, name, params)

    result = workflow.result(wait=True, poll_interval=5)
    return {
        "uuid": workflow.uuid,
        "object_data": result.data,
        "workflow_type": workflow_type,
    }


def score_candidate(cand: dict[str, str], cfg: dict[str, Any],
                    iter_n: int) -> dict[str, Any]:
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
    }

    name = f"{cfg['run_id']}_iter{iter_n:04d}_{smiles[:24]}"
    try:
        result = submit_one(smiles, cfg, name)
    except Exception as e:
        out["error"] = f"Rowan submit failed: {e}"
        out["traceback"] = traceback.format_exc()
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

    return out


def score_batch_docking(candidates: list[dict[str, str]], cfg: dict[str, Any],
                        iter_n: int) -> list[dict[str, Any]]:
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
            })
            continue

        ok, failed = check_constraints(desc, cfg.get("constraints", {}))
        scored.append({
            "smiles": smiles,
            "parent_smiles": parent,
            "design_note": note,
            "scored_at": state.now_iso(),
            "metrics": dict(desc),
            "satisfies_constraints": ok,
            "constraint_failures": failed,
            "score": None,
        })
        smiles_list.append(smiles)
        submitted_indexes.append(i)

    if not smiles_list:
        return scored

    params = _workflow_params(cfg)
    name = f"{cfg['run_id']}_iter{iter_n:04d}_batch_docking"
    try:
        workflow = rowan.submit_batch_docking_workflow(
            smiles_list=smiles_list,
            name=name,
            **params,
        )
        result = workflow.result(wait=True, poll_interval=5)
        object_data = result.data
        batch_scores = result.scores
    except Exception as e:
        error = f"Rowan batch_docking submit failed: {e}"
        tb = traceback.format_exc()
        for i in submitted_indexes:
            scored[i]["error"] = error
            scored[i]["traceback"] = tb
        return scored

    for i in submitted_indexes:
        out = scored[i]
        smiles = out["smiles"]
        val = batch_scores.get(smiles)
        out["workflow_uuid"] = workflow.uuid
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

    return scored


def rebuild_report(run_id: str) -> Path:
    """Regenerate the HTML report after an iteration is written."""
    from viz.build_report import build_report

    return build_report(run_id)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

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
def cli(run_id, rationale, candidates, decision, max_workers, dry_run):
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

    iter_n = (max(state.list_iterations(run_id)) + 1) if state.list_iterations(run_id) else 1

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
            try:
                report_path = rebuild_report(run_id)
            except Exception as e:
                click.echo(f"\nWarning: could not regenerate report.html: {e}", err=True)
                click.echo("Run `rowan-report --run %s` after fixing the report error." % run_id, err=True)
            else:
                click.echo(f"\nUpdated report: {report_path}")
                click.echo(f"Open it: file://{report_path.resolve()}")
            return

        wrapup_decision = "continue"
        if existing_iters:
            wrapup_decision = state.load_iter(run_id, existing_iters[-1]).get("decision", "continue")

        iter_record = {
            "iter": iter_n,
            "created_at": state.now_iso(),
            "rationale": rationale,
            "candidates": [],
            "iteration_best": None,
            "decision": wrapup_decision,
        }
        state.save_iter(run_id, iter_n, iter_record)
        click.echo(f"Saved wrap-up iteration {state.iter_path(run_id, iter_n)}")
        try:
            report_path = rebuild_report(run_id)
        except Exception as e:
            click.echo(f"\nWarning: could not regenerate report.html: {e}", err=True)
            click.echo("Run `rowan-report --run %s` after fixing the report error." % run_id, err=True)
        else:
            click.echo(f"\nUpdated report: {report_path}")
            click.echo(f"Open it: file://{report_path.resolve()}")
        return

    click.echo(f"[run {run_id}] iter {iter_n}: {len(parsed)} candidates")
    for c in parsed:
        click.echo(f"  - {c['smiles']}   ({c['design_note']})")

    if dry_run:
        click.echo("(dry-run; not calling Rowan)")
        return

    if not os.environ.get("ROWAN_API_KEY"):
        raise click.ClickException("ROWAN_API_KEY not set. Copy .env.example to .env and fill it in.")

    if cfg["workflow_type"] == "batch_docking":
        scored = score_batch_docking(parsed, cfg, iter_n)
        for i, c in enumerate(scored):
            metric_str = (f"{cfg['primary_metric']}={c['score']:.3f}"
                          if c.get("score") is not None else "FAILED")
            click.echo(f"  [{i+1}/{len(parsed)}] {c['smiles'][:40]:40s}  {metric_str}")
    else:
        # Submit in parallel.
        scored: list[dict[str, Any]] = [None] * len(parsed)  # type: ignore
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(score_candidate, p, cfg, iter_n): i
                       for i, p in enumerate(parsed)}
            for fut in as_completed(futures):
                i = futures[fut]
                try:
                    scored[i] = fut.result()
                except Exception as e:
                    scored[i] = {**parsed[i], "error": str(e), "score": None,
                                 "metrics": {}, "satisfies_constraints": False}
                c = scored[i]
                metric_str = (f"{cfg['primary_metric']}={c['score']:.3f}"
                              if c.get("score") is not None else "FAILED")
                click.echo(f"  [{i+1}/{len(parsed)}] {c['smiles'][:40]:40s}  {metric_str}")

    # Best across this iteration (constraint-aware, direction-aware)
    valid = [c for c in scored if c.get("score") is not None and c["satisfies_constraints"]]
    direction = cfg["objective_direction"]
    iter_best = None
    if valid:
        iter_best = (max(valid, key=lambda c: c["score"]) if direction == "maximize"
                     else min(valid, key=lambda c: c["score"]))

    iter_record = {
        "iter": iter_n,
        "created_at": state.now_iso(),
        "rationale": rationale,
        "candidates": scored,
        "iteration_best": ({"smiles": iter_best["smiles"], "score": iter_best["score"]}
                           if iter_best else None),
        "decision": decision,
    }
    state.save_iter(run_id, iter_n, iter_record)
    click.echo(f"\nSaved {state.iter_path(run_id, iter_n)}")

    best = state.best_so_far(run_id)
    if best:
        click.echo(f"Best so far: {best['smiles']}  "
                   f"{cfg['primary_metric']}={best['score']:.3f}  (iter {best['iter']})")

    try:
        report_path = rebuild_report(run_id)
    except Exception as e:
        click.echo(f"\nWarning: could not regenerate report.html: {e}", err=True)
        click.echo("Run `rowan-report --run %s` after fixing the report error." % run_id, err=True)
    else:
        click.echo(f"\nUpdated report: {report_path}")
        click.echo(f"Open it: file://{report_path.resolve()}")


if __name__ == "__main__":
    cli()
