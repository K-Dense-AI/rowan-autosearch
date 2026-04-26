"""Optional brainstorming helper. The agent is the primary proposer — this
just generates RDKit-based mutation candidates when it wants more ideas.

Strategies:
  bioisostere : swap common functional groups (COOH↔tetrazole, etc.)
  brics       : decompose with BRICS and recombine fragments
  scan-subst  : enumerate aromatic substitutions (-H to -F/-Cl/-OMe/-CN/-CF3)
"""
from __future__ import annotations

import json
import random
from typing import Iterable

import click


BIOISOSTERES: list[tuple[str, str]] = [
    # (SMARTS pattern to match, replacement SMILES fragment)
    ("[CX3](=O)[OX2H1]", "c1nnn[nH]1"),    # COOH -> tetrazole
    ("[CX3](=O)[OX2H1]", "S(=O)(=O)N"),    # COOH -> sulfonamide
    ("[CX3](=O)[NX3H2]", "S(=O)(=O)N"),    # primary amide -> sulfonamide
    ("[OX2H1]", "N"),                       # OH -> NH2
    ("[NX3H2]", "O"),                       # NH2 -> OH
    ("c1ccccc1", "c1ccncc1"),               # benzene -> pyridine
    ("c1ccccc1", "c1cncnc1"),               # benzene -> pyrimidine
    ("[CH3]", "F"),                         # methyl -> fluorine
    ("[CH3]", "C(F)(F)F"),                  # methyl -> CF3
]

AROMATIC_SUBSTS = ["F", "Cl", "Br", "C", "OC", "C(F)(F)F", "C#N", "[N+](=O)[O-]", "N"]


def bioisostere_variants(smiles: str, n: int = 8) -> list[tuple[str, str]]:
    """Return up to n (new_smiles, note) pairs."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    out: list[tuple[str, str]] = []
    seen: set[str] = {Chem.MolToSmiles(mol)}
    for patt, repl in BIOISOSTERES:
        patt_mol = Chem.MolFromSmarts(patt)
        repl_mol = Chem.MolFromSmiles(repl)
        if patt_mol is None or repl_mol is None or not mol.HasSubstructMatch(patt_mol):
            continue
        try:
            products = Chem.ReplaceSubstructs(mol, patt_mol, repl_mol, replaceAll=False)
        except Exception:
            continue
        for p in products:
            try:
                Chem.SanitizeMol(p)
            except Exception:
                continue
            s = Chem.MolToSmiles(p)
            if s in seen or "." in s:  # skip fragmented products
                continue
            seen.add(s)
            out.append((s, f"bioisostere: {patt} -> {repl}"))
            if len(out) >= n:
                return out
    return out


def aromatic_substitution_variants(smiles: str, n: int = 8) -> list[tuple[str, str]]:
    """Replace one aromatic CH with each of AROMATIC_SUBSTS."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    out: list[tuple[str, str]] = []
    seen: set[str] = {Chem.MolToSmiles(mol)}
    aromatic_h_atoms = [a.GetIdx() for a in mol.GetAtoms()
                         if a.GetIsAromatic() and a.GetTotalNumHs() > 0]
    random.shuffle(aromatic_h_atoms)
    for sub in AROMATIC_SUBSTS:
        for idx in aromatic_h_atoms:
            rwm = Chem.RWMol(mol)
            sub_mol = Chem.MolFromSmiles(sub)
            if sub_mol is None:
                continue
            offset = rwm.GetNumAtoms()
            for atom in sub_mol.GetAtoms():
                rwm.AddAtom(atom)
            for b in sub_mol.GetBonds():
                rwm.AddBond(b.GetBeginAtomIdx() + offset,
                            b.GetEndAtomIdx() + offset, b.GetBondType())
            rwm.AddBond(idx, offset, Chem.BondType.SINGLE)
            try:
                m2 = rwm.GetMol()
                Chem.SanitizeMol(m2)
                s = Chem.MolToSmiles(m2)
            except Exception:
                continue
            if s in seen:
                continue
            seen.add(s)
            out.append((s, f"add -{sub} on aromatic ring"))
            if len(out) >= n:
                return out
            break  # one position per substituent
    return out


def brics_variants(smiles: str, n: int = 8) -> list[tuple[str, str]]:
    from rdkit import Chem
    from rdkit.Chem import BRICS

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    frags = list(BRICS.BRICSDecompose(mol))
    if len(frags) < 2:
        return []
    builder = BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in frags if Chem.MolFromSmiles(f)])
    out: list[tuple[str, str]] = []
    seen: set[str] = {Chem.MolToSmiles(mol)}
    for i, p in enumerate(builder):
        if i > 200:
            break
        try:
            Chem.SanitizeMol(p)
        except Exception:
            continue
        s = Chem.MolToSmiles(p)
        if s in seen:
            continue
        seen.add(s)
        out.append((s, "BRICS recombination"))
        if len(out) >= n:
            return out
    return out


@click.command()
@click.option("--smiles", required=True, help="Parent SMILES.")
@click.option("--strategy",
              type=click.Choice(["bioisostere", "scan-subst", "brics", "all"]),
              default="all")
@click.option("--n", default=8, type=int, help="Max variants per strategy.")
@click.option("--seed", default=0, type=int)
def cli(smiles, strategy, n, seed):
    """Print candidate SMILES variants for the given parent.

    Output format is one '<smiles>|<parent>|<note>' per line — pipe directly
    into `rowan-score --candidate ...` after the agent has reviewed them.
    """
    random.seed(seed)
    fns = {
        "bioisostere": bioisostere_variants,
        "scan-subst": aromatic_substitution_variants,
        "brics": brics_variants,
    }
    if strategy == "all":
        keys = list(fns)
    else:
        keys = [strategy]

    variants: list[tuple[str, str]] = []
    for k in keys:
        variants.extend(fns[k](smiles, n=n))

    seen: set[str] = set()
    for s, note in variants:
        if s in seen:
            continue
        seen.add(s)
        click.echo(f"{s}|{smiles}|{note}")


if __name__ == "__main__":
    cli()
