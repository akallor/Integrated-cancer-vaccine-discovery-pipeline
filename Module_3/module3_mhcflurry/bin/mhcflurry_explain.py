#!/usr/bin/env python3
"""
mhcflurry_explain.py  --  Module 3 explainability layer for the neoantigen pipeline.

Adds two interpretability outputs on top of MHCflurry HLA-I predictions:

  (a) POSITIONAL SHAP  -- per-residue attribution showing which positions of a
      peptide MHCflurry relies on to decide strong- vs weak-binding, computed with
      a KernelExplainer over the MHCflurry affinity predictor treated as a black box.
      Positional attributions are compared against the *observed* binding motif
      (empirical PWM built from the sample's own strong binders for that allele/length),
      i.e. "does the model weight the same anchor residues the real HLA groove selects?".

  (b) ISM "WEIGHTS MATRIX" -- in-silico saturation mutagenesis. For each selected
      peptide we mutate every position to all 20 amino acids and record the change in
      the model's binding score, giving an L x 20 matrix of the effective per-position
      weights the model assigned. This is the model-internal motif and is plotted as a
      heatmap and as a sequence logo alongside the empirical (real) motif.

Peptide selection: 2 peptides from each of the four strata
      {group} x {binder call}  ->  CANON/SB, CANON/WB, NONCANON/SB, NONCANON/WB  (8 total),
sampled reproducibly with --seed.

Backends:
  --backend mhcflurry  : uses the real Class1AffinityPredictor (run this on EC2).
  --backend surrogate  : a PWM-based stand-in with the identical call signature, used
                         only to validate the plumbing where MHCflurry is unavailable.

The two backends are interchangeable because everything downstream only ever calls
  predictor.predict_affinity(peptides: list[str], allele: str) -> np.ndarray[nM]

Author: (pipeline Module 3)  --  tested end-to-end in surrogate mode.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

AA = "ACDEFGHIKLMNPQRSTVWY"           # 20 canonical amino acids, fixed order
AA_INDEX = {a: i for i, a in enumerate(AA)}


# --------------------------------------------------------------------------------------
# Affinity <-> score transform (standard NetMHC log50k), higher score == stronger binder
# --------------------------------------------------------------------------------------
def aff_to_score(aff_nm: np.ndarray) -> np.ndarray:
    """1 - log50k(affinity). ~1.0 for a 1 nM binder, ~0.0 for >=50000 nM. Clipped to [0,1]."""
    aff_nm = np.clip(np.asarray(aff_nm, dtype=float), 1e-3, 5e4)
    return np.clip(1.0 - np.log(aff_nm) / np.log(5e4), 0.0, 1.0)


# --------------------------------------------------------------------------------------
# Backends -- both expose predict_affinity(peptides, allele) -> np.ndarray of nM
# --------------------------------------------------------------------------------------
class MHCflurryBackend:
    """Wraps the real MHCflurry 2.0 Class1AffinityPredictor. Use on EC2."""

    def __init__(self):
        from mhcflurry import Class1PresentationPredictor  # imported lazily
        # The presentation predictor bundles the affinity predictor we want to explain.
        self._pres = Class1PresentationPredictor.load()
        self._aff = self._pres.affinity_predictor
        self.supported_alleles = set(self._aff.supported_alleles)

    def predict_affinity(self, peptides, allele):
        # Class1AffinityPredictor.predict returns nM affinities for a single allele.
        return np.asarray(
            self._aff.predict(peptides=list(peptides), allele=allele, throw=False),
            dtype=float,
        )


class SurrogateBackend:
    """
    PWM-based stand-in with the SAME signature as MHCflurryBackend, so the entire
    SHAP / ISM / motif / plotting stack can be validated without TensorFlow or the
    MHCflurry model weights. Anchor preferences are allele-specific and biologically
    shaped (P2 and P-Omega dominant for HLA-I) so the outputs look realistic.

    NOT a predictor of real affinities -- validation only.
    """

    # rough, illustrative anchor preferences per allele (position: {aa: bonus})
    _MOTIFS = {
        "HLA-B*15:01": {2: {"Q": 3.0, "K": 2.0, "L": 1.0}, -1: {"Y": 3.0, "F": 2.5, "W": 1.5}},
        "HLA-C*12:03": {2: {"A": 2.0, "L": 2.0, "F": 1.5}, -1: {"L": 3.0, "F": 2.0, "Y": 1.5}},
        "HLA-A*02:01": {2: {"L": 3.0, "M": 2.5, "I": 1.5}, -1: {"V": 3.0, "L": 2.5, "A": 1.5}},
        "HLA-A*01:01": {3: {"D": 3.0, "E": 2.0}, -1: {"Y": 3.0, "F": 2.0}},
    }
    _DEFAULT = {2: {"L": 2.0, "I": 1.5}, -1: {"Y": 2.0, "L": 1.5, "F": 1.5}}

    def __init__(self, seed: int = 0):
        self._rng = np.random.default_rng(seed)
        # small random per-(pos,aa) background so non-anchor positions still matter a little
        self._bg = {}

    def _bg_matrix(self, allele, L):
        key = (allele, L)
        if key not in self._bg:
            self._bg[key] = self._rng.normal(0, 0.15, size=(L, 20))
        return self._bg[key]

    def predict_affinity(self, peptides, allele):
        motif = self._MOTIFS.get(allele, self._DEFAULT)
        out = []
        for pep in peptides:
            L = len(pep)
            bg = self._bg_matrix(allele, L)
            s = 0.0
            for i, aa in enumerate(pep):
                if aa in AA_INDEX:
                    s += bg[i, AA_INDEX[aa]]
            for pos, prefs in motif.items():
                idx = pos if pos >= 0 else L + pos
                if 0 <= idx < L:
                    s += prefs.get(pep[idx], 0.0)
            # map summed score -> pseudo-affinity nM (higher s -> lower nM)
            aff = float(50000.0 * np.exp(-1.15 * s))
            out.append(np.clip(aff, 1.0, 50000.0))
        return np.asarray(out, dtype=float)


def load_backend(name: str):
    if name == "mhcflurry":
        return MHCflurryBackend()
    if name == "surrogate":
        return SurrogateBackend()
    raise ValueError(f"unknown backend {name!r}")


# --------------------------------------------------------------------------------------
# Peptide selection: 2 per (group x binder-call) stratum
# --------------------------------------------------------------------------------------
@dataclass
class Selected:
    peptide: str
    allele: str
    group: str          # CANON / NONCANON
    call: str           # SB / WB
    affinity: float
    row: pd.Series


def _norm_group(v: str) -> str:
    s = str(v).strip().lower()
    if s.startswith("non") or "noncanon" in s:
        return "NONCANON"
    return "CANON"


def select_peptides(df, group_col, allele_col, call_col, aff_col,
                    n_per_stratum=2, seed=0, calls=("SB", "WB")):
    rng = np.random.default_rng(seed)
    g = df.copy()
    g["_group"] = g[group_col].map(_norm_group)
    g["_call"] = g[call_col].astype(str).str.upper().str.strip()

    chosen = []
    for grp in ("CANON", "NONCANON"):
        for call in calls:
            cell = g[(g["_group"] == grp) & (g["_call"] == call)]
            # de-duplicate on (peptide, allele) so we don't pick the same pMHC twice
            cell = cell.drop_duplicates(subset=[ "peptide", allele_col])
            if len(cell) == 0:
                print(f"  [warn] no rows for stratum {grp}/{call} -- skipping", file=sys.stderr)
                continue
            take = min(n_per_stratum, len(cell))
            idx = rng.choice(cell.index.to_numpy(), size=take, replace=False)
            for i in idx:
                r = cell.loc[i]
                chosen.append(Selected(
                    peptide=str(r["peptide"]).strip().upper(),
                    allele=str(r[allele_col]).strip(),
                    group=grp, call=call,
                    affinity=float(r[aff_col]),
                    row=r,
                ))
    return chosen


# --------------------------------------------------------------------------------------
# Background peptide pool for SHAP masking (per allele / length)
# --------------------------------------------------------------------------------------
def build_background_pool(df, allele_col, allele, L, n=120, seed=0,
                          mode="dataset", call_col=None):
    """
    Returns an (n, L) integer matrix of AA indices used by KernelExplainer to
    substitute 'off' positions. mode='dataset' samples real peptides of the same
    allele+length from the sample's own predictions (SHAP then explains a peptide
    relative to the natural repertoire for that allele). mode='random' uses a
    uniform amino-acid background (explains relative to random sequence).
    """
    rng = np.random.default_rng(seed)
    peps = []
    if mode == "dataset":
        sub = df[(df[allele_col].astype(str).str.strip() == allele)]
        sub = sub[sub["peptide"].astype(str).str.len() == L]
        peps = sub["peptide"].astype(str).str.upper().tolist()
        # keep only clean L-mers over the 20 canonical AAs
        peps = [p for p in peps if len(p) == L and set(p) <= set(AA)]
    if len(peps) < 20:  # fall back to random if the pool is too thin
        peps = ["".join(rng.choice(list(AA), size=L)) for _ in range(max(n, 60))]
    sel = rng.choice(len(peps), size=min(n, len(peps)), replace=len(peps) < n)
    mat = np.array([[AA_INDEX[a] for a in peps[k]] for k in sel], dtype=int)
    return mat


def _decode(mat_row):
    return "".join(AA[int(i)] for i in mat_row)


# --------------------------------------------------------------------------------------
# (a) Positional SHAP
# --------------------------------------------------------------------------------------
def positional_shap(backend, peptide, allele, background_mat, nsamples="auto"):
    """
    KernelExplainer over an integer-encoded peptide (one feature per position).
    Perturbing a feature == swapping that residue for the corresponding residue of a
    background peptide. Returns per-position SHAP on the binding SCORE (higher score =
    stronger binder), so positive SHAP == residue pushes toward strong binding.
    """
    import shap

    L = len(peptide)
    x = np.array([[AA_INDEX[a] for a in peptide]], dtype=float)   # (1, L)
    bg = background_mat.astype(float)

    def f(mat):
        peps = [_decode(np.rint(row).astype(int)) for row in np.atleast_2d(mat)]
        return aff_to_score(backend.predict_affinity(peps, allele))

    explainer = shap.KernelExplainer(f, bg, silent=True)
    sv = explainer.shap_values(x, nsamples=nsamples, silent=True)
    sv = np.array(sv).reshape(-1)[:L]
    base = float(np.mean(aff_to_score(backend.predict_affinity(
        [_decode(r) for r in np.rint(bg).astype(int)], allele))))
    wt_score = float(aff_to_score(backend.predict_affinity([peptide], allele))[0])
    return sv, base, wt_score


# --------------------------------------------------------------------------------------
# (b) In-silico saturation mutagenesis -> the model "weights matrix"
# --------------------------------------------------------------------------------------
def ism_matrix(backend, peptide, allele):
    """
    L x 20 matrix. Entry (i, a) = binding score of the peptide with position i mutated
    to amino acid a. The wild-type column is the diagonal. Delta version (mut - wt) is
    returned too -- that is the effective per-position weight the model assigns.
    """
    L = len(peptide)
    variants, coords = [], []
    for i in range(L):
        for a_i, a in enumerate(AA):
            v = peptide[:i] + a + peptide[i + 1:]
            variants.append(v)
            coords.append((i, a_i))
    scores = aff_to_score(backend.predict_affinity(variants, allele))
    M = np.zeros((L, 20), dtype=float)
    for (i, a_i), s in zip(coords, scores):
        M[i, a_i] = s
    wt = aff_to_score(backend.predict_affinity([peptide], allele))[0]
    delta = M - wt
    return M, delta, float(wt)


# --------------------------------------------------------------------------------------
# Empirical ("real") motif: PWM from the sample's own strong binders (allele + length)
# --------------------------------------------------------------------------------------
def empirical_pwm(df, allele_col, call_col, allele, L, strong_calls=("SB",),
                  pseudocount=0.5):
    sub = df[(df[allele_col].astype(str).str.strip() == allele)
             & (df[call_col].astype(str).str.upper().str.strip().isin(strong_calls))]
    peps = [p.upper() for p in sub["peptide"].astype(str) if len(p) == L and set(p.upper()) <= set(AA)]
    counts = np.full((L, 20), pseudocount, dtype=float)
    for p in peps:
        for i, a in enumerate(p):
            counts[i, AA_INDEX[a]] += 1.0
    pwm = counts / counts.sum(axis=1, keepdims=True)
    return pwm, len(peps)


def pwm_to_bits(pwm):
    """Convert a probability PWM to information-content (bits) heights for a logo."""
    L = pwm.shape[0]
    ic = np.log2(20) + np.sum(np.where(pwm > 0, pwm * np.log2(pwm + 1e-12), 0.0), axis=1)
    ic = np.clip(ic, 0, np.log2(20))
    return pwm * ic[:, None]


def ism_to_pwm(M, temp=6.0):
    """Turn a positive ISM score matrix into a softmax preference PWM per position."""
    Z = M - M.max(axis=1, keepdims=True)
    P = np.exp(temp * Z)
    return P / P.sum(axis=1, keepdims=True)


def model_motif_by_sampling(backend, allele, L, n_random=3000, top_frac=0.02,
                            seed=0, pseudocount=0.5):
    """
    The model's *intrinsic* binding motif, computed the same way the observed motif is:
    as residue frequencies among strong binders -- except here the "strong binders" are
    the peptides the MODEL scores highest. Score n_random random L-mers, keep the top
    fraction, and build a PWM. Directly comparable to empirical_pwm() and independent of
    the particular selected peptides.
    """
    rng = np.random.default_rng(seed)
    peps = ["".join(rng.choice(list(AA), size=L)) for _ in range(n_random)]
    scores = aff_to_score(backend.predict_affinity(peps, allele))
    k = max(20, int(n_random * top_frac))
    top = np.argsort(scores)[::-1][:k]
    counts = np.full((L, 20), pseudocount, dtype=float)
    for j in top:
        for i, a in enumerate(peps[j]):
            counts[i, AA_INDEX[a]] += 1.0
    pwm = counts / counts.sum(axis=1, keepdims=True)
    return pwm, k
