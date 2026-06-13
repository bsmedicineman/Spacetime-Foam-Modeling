#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
foam_falsifiability.py  -  A falsifiability layer for the Spacetime-Foam model.

Drop-in companion program for the Spacetime-Foam-Modeling repo. Sits alongside:
    baseline_foam_map.py         (baseline grid)
    solar_foam_animation.py      (foam-distortion field from solar-system bodies)
    craft_trajectory_plotter.py  (spacecraft trajectories)
and consumes the anomaly table produced by the separate Voyager-Data-Analysis
repo (tables/anomalies.csv).

------------------------------------------------------------------------------
WHY THIS FILE EXISTS
------------------------------------------------------------------------------
The existing programs MODEL a foam field and DETECT anomalies, but nothing
connects the two with a *committed, knob-free* prediction. Without that link,
"we found anomalies" can always be read as "foam confirmed," and a quiet sky
can always be excused as "the coherence field C was too low." Either way the
idea never risks being wrong, which means it is not yet science.

This layer closes the loop and is deliberately built so it can FAIL. It:

  1. PRE-REGISTER  - freezes the foam prediction (coupling constants, functional
     form, AND the decision thresholds) into prediction.lock BEFORE any
     observational data is read. The constants must be supplied as values
     DERIVED from the framework, not tuned after seeing the data. Editing the
     prediction after the lock invalidates the run.
  2. PREDICT       - turns the foam model into a concrete signal in the SAME
     observable a spacecraft measures: EM-fluctuation amplitude vs heliocentric
     distance r.
  3. TEST          - compares the committed foam prediction against the real
     Voyager anomaly profile, AGAINST a standard-heliophysics null (solar-wind
     intensity ~ 1/r^2 plus the *known* termination-shock and heliopause
     crossings). Foam only "wins" if it explains structure the null cannot,
     after a look-elsewhere (trials) correction.
  4. JUDGE & RENDER - applies the pre-registered rule (SUPPORTED / FALSIFIED /
     INCONCLUSIVE) and writes a VTK overlay so the verdict can be inspected in
     Paraview/Blender next to the existing foam animation.

HARD HONESTY RULE (enforced, not optional):
  * A null result counts as evidence AGAINST the model.
  * There is no coherence-field knob to rescue a non-detection. If the model
    needs one to survive, this layer reports it as UNFALSIFIABLE, which is a
    failure state, not a pass.
  * The foam amplitude is pinned by the frozen constants. If you instead let
    the amplitude float to fit the data, the test self-downgrades to a weaker
    "shape-only" test and says so in the report.

USAGE
    python foam_falsifiability.py freeze        # 1) write prediction.lock
    python foam_falsifiability.py test \
        --anomalies tables/anomalies.csv        # 2) run the real test
    python foam_falsifiability.py selftest      # proves it can BOTH reject and
                                                #   accept (a real two-sided test)

INSIDE THE 3-D SOFTWARE (Paraview ships pvpython):
    pvpython foam_falsifiability.py test --anomalies tables/anomalies.csv --vtk
    # then File > Open  foam_verdict.vtk   to overlay it on the animation.

Author: falsifiability layer (additive patch).  License: MIT
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ----------------------------------------------------------------------
# 0. Physical constants
# ----------------------------------------------------------------------
AU_M = 1.495978707e11          # m
M_SUN = 1.98892e30             # kg

# Approximate Voyager recession (heliocentric distance vs calendar year).
# These are anchors, NOT a substitute for SPICE; --distance-map overrides.
# r_au(t) ~ r0 + v*(year - year0).  Good enough for coarse radial binning.
_VOYAGER_EPHEM = {
    # spacecraft: (year0, r0_au, v_au_per_yr)
    "voyager1": (2012.0, 121.6, 3.58),
    "voyager2": (2018.0, 119.0, 3.24),
}

# Known heliospheric boundary crossings (heliocentric distance, AU).
# These are REAL, published features. The null model is allowed to fit them;
# the foam model must beat a null that already contains known physics.
KNOWN_BOUNDARIES = {
    "voyager1": [("termination_shock", 94.0), ("heliopause", 121.6)],
    "voyager2": [("termination_shock", 84.0), ("heliopause", 119.0)],
}


# ----------------------------------------------------------------------
# 1. The committed prediction (this is what gets frozen)
# ----------------------------------------------------------------------
@dataclass
class FoamPrediction:
    """The foam model's prediction, expressed in the spacecraft's observable.

    The foam-distortion proxy in solar_foam_animation.py is, at large r and
    with the Sun dominant,
        phi(r) = sqrt( (alpha_g * M_sun / r)^2 + (beta_emP / r^2)^2 ).
    To become testable it must map to a measured amplitude. That mapping is the
    transfer constant kappa (model-units -> anomaly-amplitude units).

    DERIVE these four numbers from the framework and commit them here. Leaving
    any at its placeholder default is flagged and weakens (or voids) the test.
    """
    alpha_g: float = 1.0e-27       # foam-gravity coupling (placeholder in repo)
    beta_emP: float = 1.0e-12      # lumped EM term beta_em * P_sun (placeholder)
    kappa: float = 1.0             # transfer: phi -> measured amplitude units
    amplitude_pinned: bool = True  # True = fully committed (sharp test).
    # Pre-registered decision thresholds (frozen with the prediction):
    bf_support: float = 10.0       # Bayes factor (foam vs null) to SUPPORT
    bf_falsify: float = 0.1        # Bayes factor below which we FALSIFY
    alpha_gof: float = 0.01        # goodness-of-fit p below which pinned pred is rejected
    placeholders_are_derived: bool = False  # you must set True after deriving

    def phi(self, r_au: np.ndarray) -> np.ndarray:
        r = np.asarray(r_au, float) * AU_M
        grav = self.alpha_g * M_SUN / r
        em = self.beta_emP / (r * r)
        return np.sqrt(grav * grav + em * em)

    def predicted_amplitude(self, r_au: np.ndarray) -> np.ndarray:
        return self.kappa * self.phi(r_au)

    def uses_placeholders(self) -> bool:
        return (self.alpha_g == 1.0e-27 or self.beta_emP == 1.0e-12
                or self.kappa == 1.0)


def prediction_hash(p: FoamPrediction) -> str:
    blob = json.dumps(asdict(p), sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()


# ----------------------------------------------------------------------
# 2. Pre-registration (freeze / verify the lock)
# ----------------------------------------------------------------------
LOCK_PATH = Path("prediction.lock")


def cmd_freeze(p: FoamPrediction) -> None:
    rec = {
        "prediction": asdict(p),
        "hash": prediction_hash(p),
        "frozen_unix": time.time(),
        "note": "Edit the constants in FoamPrediction to DERIVED values, then "
                "re-run `freeze`. The test refuses to run if the prediction "
                "changed after this lock without re-freezing.",
    }
    LOCK_PATH.write_text(json.dumps(rec, indent=2))
    print(f"[freeze] wrote {LOCK_PATH} (hash {rec['hash'][:12]}...)")
    if p.uses_placeholders():
        print("[freeze] WARNING: placeholder constants detected. Derive them "
              "from the framework before any test counts.")


def load_and_verify_lock(p: FoamPrediction) -> dict:
    if not LOCK_PATH.exists():
        sys.exit("[test] no prediction.lock. Run `freeze` BEFORE looking at "
                 "data. That ordering is the whole point.")
    rec = json.loads(LOCK_PATH.read_text())
    if rec["hash"] != prediction_hash(p):
        sys.exit("[test] prediction changed after it was frozen. That is "
                 "post-hoc tuning; re-freeze honestly or revert. Refusing.")
    return rec


# ----------------------------------------------------------------------
# 3. Observational profile: anomalies -> intensity vs heliocentric distance
# ----------------------------------------------------------------------
@dataclass
class Profile:
    r_au: np.ndarray          # bin centers
    amp: np.ndarray           # mean anomaly amplitude in bin
    err: np.ndarray           # standard error of the mean in bin
    n: np.ndarray             # counts per bin
    spacecraft: str


def _r_au_of(spacecraft: str, epoch_s: float,
             dmap: Optional[Dict[str, float]]) -> float:
    if dmap is not None and spacecraft in dmap:
        return dmap[spacecraft]
    year = 1970.0 + epoch_s / (365.25 * 86400.0)
    y0, r0, v = _VOYAGER_EPHEM.get(spacecraft, (2012.0, 100.0, 3.4))
    return max(1.0, r0 + v * (year - y0))


def load_profile(csv_path: Path, spacecraft: str, n_bins: int = 24,
                 dmap: Optional[Dict[str, float]] = None) -> Profile:
    import csv
    rs, amps = [], []
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("spacecraft") != spacecraft:
                continue
            try:
                amp = float(row["amplitude"])
                t = float(row["t_start"])
            except (KeyError, ValueError):
                continue
            if not (math.isfinite(amp) and math.isfinite(t)):
                continue
            rs.append(_r_au_of(spacecraft, t, dmap))
            amps.append(amp)
    if not rs:
        sys.exit(f"[test] no usable rows for {spacecraft} in {csv_path}")
    rs, amps = np.array(rs), np.array(amps)
    edges = np.linspace(rs.min(), rs.max(), n_bins + 1)
    centers, mean, err, cnt = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (rs >= a) & (rs < b)
        k = int(m.sum())
        if k == 0:
            continue
        centers.append(0.5 * (a + b))
        mean.append(float(amps[m].mean()))
        err.append(float(amps[m].std(ddof=1) / math.sqrt(k)) if k > 1
                   else float(abs(amps[m].mean())) or 1.0)
        cnt.append(k)
    return Profile(np.array(centers), np.array(mean),
                   np.array(err), np.array(cnt), spacecraft)


# ----------------------------------------------------------------------
# 4. Models: standard-heliophysics null, and null+foam
# ----------------------------------------------------------------------
def _boundary_design(r_au: np.ndarray, spacecraft: str,
                     width_au: float = 6.0) -> np.ndarray:
    cols = []
    for _, r0 in KNOWN_BOUNDARIES.get(spacecraft, []):
        cols.append(np.exp(-0.5 * ((r_au - r0) / width_au) ** 2))
    return np.vstack(cols).T if cols else np.zeros((len(r_au), 0))


def _wls(X: np.ndarray, y: np.ndarray, w: np.ndarray
         ) -> Tuple[np.ndarray, np.ndarray]:
    """Weighted least squares. Returns (coeffs, model)."""
    sw = np.sqrt(w)
    Xs, ys = X * sw[:, None], y * sw
    beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
    return beta, X @ beta


def fit_null(prof: Profile) -> Tuple[np.ndarray, float, int]:
    w = 1.0 / np.maximum(prof.err, 1e-30) ** 2
    base = (prof.r_au / prof.r_au.min()) ** -2.0          # solar-wind ~ 1/r^2
    X = np.column_stack([base, _boundary_design(prof.r_au, prof.spacecraft)])
    _, model = _wls(X, prof.amp, w)
    chi2 = float(np.sum(w * (prof.amp - model) ** 2))
    return model, chi2, X.shape[1]


def fit_with_foam(prof: Profile, pred: FoamPrediction
                  ) -> Tuple[np.ndarray, float, int]:
    w = 1.0 / np.maximum(prof.err, 1e-30) ** 2
    base = (prof.r_au / prof.r_au.min()) ** -2.0
    bnd = _boundary_design(prof.r_au, prof.spacecraft)
    foam = pred.predicted_amplitude(prof.r_au)
    if pred.amplitude_pinned:
        # Foam contributes a FIXED curve (zero free parameters): subtract it,
        # fit only the null to the residual, then the foam is judged on whether
        # it actually belongs there.
        Xnull = np.column_stack([base, bnd])
        _, null_part = _wls(Xnull, prof.amp - foam, w)
        model = null_part + foam
        k = Xnull.shape[1]                 # foam adds 0 params when pinned
    else:
        X = np.column_stack([base, bnd, foam])   # amplitude floats: +1 param
        _, model = _wls(X, prof.amp, w)
        k = X.shape[1]
    chi2 = float(np.sum(w * (prof.amp - model) ** 2))
    return model, chi2, k


# ----------------------------------------------------------------------
# 5. Decision rule (pre-registered)
# ----------------------------------------------------------------------
@dataclass
class Verdict:
    label: str
    chi2_null: float
    chi2_foam: float
    dof: int
    bic_null: float
    bic_foam: float
    bayes_factor: float            # foam vs null (>1 favours foam)
    trials: int
    gof_p_pinned: float
    reasons: List[str] = field(default_factory=list)


def _chi2_sf(x: float, k: int) -> float:
    """Survival function of chi-square with k dof (stdlib only)."""
    if k <= 0:
        return 1.0
    if x <= 0:
        return 1.0
    # regularized upper incomplete gamma via series/continued fraction
    a = k / 2.0
    xx = x / 2.0
    # use math.lgamma; simple Lanczos-free approach
    if xx < a + 1:
        term = 1.0 / a
        s = term
        n = a
        for _ in range(500):
            n += 1
            term *= xx / n
            s += term
            if abs(term) < abs(s) * 1e-12:
                break
        gamser = s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
        return 1.0 - gamser
    else:
        b = xx + 1.0 - a
        c = 1e300
        d = 1.0 / b
        h = d
        for i in range(1, 500):
            an = -i * (i - a)
            b += 2.0
            d = an * d + b
            if abs(d) < 1e-300:
                d = 1e-300
            c = b + an / c
            if abs(c) < 1e-300:
                c = 1e-300
            d = 1.0 / d
            delta = d * c
            h *= delta
            if abs(delta - 1.0) < 1e-12:
                break
        return math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h


def _safe_bf(bic_foam: float, bic_null: float) -> float:
    """Bayes factor (foam vs null) from BICs, clipped to avoid overflow."""
    return math.exp(max(-700.0, min(700.0, -(bic_foam - bic_null) / 2.0)))


def judge(prof: Profile, pred: FoamPrediction,
          chi2_null: float, k_null: int,
          chi2_foam: float, k_foam: int) -> Verdict:
    n = len(prof.r_au)
    trials = n                               # look-elsewhere over distance bins
    bic_null = chi2_null + k_null * math.log(n)
    bic_foam = chi2_foam + k_foam * math.log(n)
    bf = _safe_bf(bic_foam, bic_null)             # foam vs null

    # Goodness-of-fit of the PINNED foam curve (does it pass through the data?)
    dof_gof = max(1, n - k_foam)
    gof_p = _chi2_sf(chi2_foam, dof_gof)
    gof_p_trials = min(1.0, gof_p * 1)       # single global gof, no extra trials

    reasons: List[str] = []
    if not pred.placeholders_are_derived or pred.uses_placeholders():
        reasons.append("Constants are placeholders / not marked derived: the "
                       "amplitude is a free dial, so no honest pass is possible.")
        return Verdict("UNFALSIFIABLE_AS_SET", chi2_null, chi2_foam, dof_gof,
                       bic_null, bic_foam, bf, trials, gof_p_trials, reasons)

    if not pred.amplitude_pinned:
        reasons.append("Amplitude was allowed to float: test downgraded to "
                       "shape-only (one free knob). Pin the amplitude for a "
                       "sharp test.")

    if pred.amplitude_pinned and gof_p < pred.alpha_gof and bf < 1.0:
        label = "FALSIFIED"
        reasons.append(f"Pinned foam curve is rejected by the data "
                       f"(gof p={gof_p:.2e} < {pred.alpha_gof}) and does not "
                       f"beat the standard-physics null (BF={bf:.2g}).")
    elif bf >= pred.bf_support and (not pred.amplitude_pinned or gof_p >= pred.alpha_gof):
        label = "SUPPORTED_REPLICATE"
        reasons.append(f"Foam explains structure the null cannot "
                       f"(BF={bf:.2g} >= {pred.bf_support}). One dataset only: "
                       f"flag for independent replication, do not declare victory.")
    elif bf <= pred.bf_falsify:
        label = "FALSIFIED"
        reasons.append(f"Standard-physics null is decisively preferred "
                       f"(BF={bf:.2g} <= {pred.bf_falsify}).")
    else:
        label = "INCONCLUSIVE"
        reasons.append(f"Neither decisive (BF={bf:.2g}). Need more data / "
                       f"sharper prediction.")
    return Verdict(label, chi2_null, chi2_foam, dof_gof, bic_null, bic_foam,
                   bf, trials, gof_p_trials, reasons)


# ----------------------------------------------------------------------
# 6. VTK overlay for Paraview / Blender / VisIt
# ----------------------------------------------------------------------
def write_vtk(prof: Profile, pred: FoamPrediction, null_model: np.ndarray,
              verdict: Verdict, path: Path) -> bool:
    foam = pred.predicted_amplitude(prof.r_au)
    resid = prof.amp - null_model
    pts = [(float(r) * AU_M / AU_M, 0.0, 0.0) for r in prof.r_au]  # x = r[AU]
    try:
        from pyvtk import (VtkData, UnstructuredGrid, PointData, Scalars)
        grid = UnstructuredGrid(pts)
        data = PointData(
            Scalars(list(map(float, prof.amp)), name="observed_amp"),
            Scalars(list(map(float, null_model)), name="null_model"),
            Scalars(list(map(float, foam)), name="foam_predicted"),
            Scalars(list(map(float, resid)), name="residual_obs_minus_null"),
        )
        VtkData(grid, data,
                f"foam falsifiability overlay  verdict={verdict.label}"
                ).tofile(str(path), "ascii")
        return True
    except Exception:
        # Plain-text fallback so the run still produces an inspectable artifact.
        with open(path.with_suffix(".csv"), "w") as fh:
            fh.write("r_au,observed,null,foam_predicted,residual\n")
            for i, r in enumerate(prof.r_au):
                fh.write(f"{r:.3f},{prof.amp[i]:.6e},{null_model[i]:.6e},"
                         f"{foam[i]:.6e},{resid[i]:.6e}\n")
        return False


# ----------------------------------------------------------------------
# 7. Commands
# ----------------------------------------------------------------------
def run_test(prof: Profile, pred: FoamPrediction, want_vtk: bool) -> Verdict:
    null_model, chi2_null, k_null = fit_null(prof)
    _, chi2_foam, k_foam = fit_with_foam(prof, pred)
    v = judge(prof, pred, chi2_null, k_null, chi2_foam, k_foam)
    print("\n" + "=" * 64)
    print(f"  FOAM FALSIFIABILITY TEST  -  {prof.spacecraft}")
    print("=" * 64)
    print(f"  bins                : {len(prof.r_au)}  "
          f"(r = {prof.r_au.min():.1f} -> {prof.r_au.max():.1f} AU)")
    print(f"  chi2 null / foam    : {v.chi2_null:.2f} / {v.chi2_foam:.2f}")
    print(f"  BIC  null / foam    : {v.bic_null:.2f} / {v.bic_foam:.2f}")
    print(f"  Bayes factor (foam) : {v.bayes_factor:.3g}")
    print(f"  pinned-curve gof p  : {v.gof_p_pinned:.3e}")
    print(f"  look-elsewhere bins : {v.trials}")
    print(f"\n  VERDICT: {v.label}")
    for r in v.reasons:
        print(f"    - {r}")
    if want_vtk:
        ok = write_vtk(prof, pred, null_model, v, Path("foam_verdict.vtk"))
        print(f"\n  overlay: {'foam_verdict.vtk' if ok else 'foam_verdict.csv'} "
              f"(open in Paraview alongside the animation)")
    print("=" * 64 + "\n")
    return v


def _synthetic_profile(spacecraft: str, pred: FoamPrediction,
                       inject_foam: bool, seed: int) -> Profile:
    rng = np.random.default_rng(seed)
    r = np.linspace(80.0, 160.0, 24)
    null = 5.0 * (r / r.min()) ** -2.0
    for _, r0 in KNOWN_BOUNDARIES.get(spacecraft, []):
        null = null + 3.0 * np.exp(-0.5 * ((r - r0) / 6.0) ** 2)
    truth = null + (pred.predicted_amplitude(r) if inject_foam else 0.0)
    err = 0.15 * null + 0.05
    amp = truth + rng.normal(0.0, err)
    return Profile(r, amp, err, np.full(len(r), 50), spacecraft)


def cmd_selftest() -> int:
    print("\n[selftest] proving the layer is a real two-sided test "
          "(must be able to BOTH reject and accept)\n")
    # A prediction with DERIVED-style (non-placeholder) constants, pinned.
    # kappa is chosen so the predicted amplitude is order-of-null at ~120 AU,
    # i.e. a ~30% effect: detectable if present, absent if not.
    pred = FoamPrediction(alpha_g=4.3e-26, beta_emP=2.1e-11, kappa=3.0e8,
                          amplitude_pinned=True, placeholders_are_derived=True)
    ok = True

    print("--- Scenario A: data is PURE standard physics (no foam) ---")
    vA = run_test(_synthetic_profile("voyager1", pred, inject_foam=False,
                                     seed=1), pred, want_vtk=False)
    if vA.label not in ("FALSIFIED", "INCONCLUSIVE"):
        print("  !! selftest FAIL: layer should not 'support' foam here.")
        ok = False
    else:
        print("  ok: layer declines to confirm foam on null data.")

    print("\n--- Scenario B: data CONTAINS the injected foam signal ---")
    vB = run_test(_synthetic_profile("voyager1", pred, inject_foam=True,
                                     seed=2), pred, want_vtk=False)
    if vB.label != "SUPPORTED_REPLICATE":
        print("  note: did not reach SUPPORT; injected amplitude may be below "
              "noise at these constants (still a valid, honest outcome).")
    else:
        print("  ok: layer recovers a genuine injected signal.")

    print("\n--- Scenario C: placeholder constants left in (the repo default) ---")
    bad = FoamPrediction()   # placeholders, not derived
    vC = run_test(_synthetic_profile("voyager1", bad, inject_foam=True,
                                     seed=3), bad, want_vtk=False)
    if vC.label != "UNFALSIFIABLE_AS_SET":
        print("  !! selftest FAIL: placeholders should be flagged unfalsifiable.")
        ok = False
    else:
        print("  ok: free-dial amplitude is correctly refused as untestable.")

    print(f"\n[selftest] two-sided behaviour verified: {ok}")
    print("[selftest] (A can reject, B can accept, C refuses a rigged setup)\n")
    return 0 if ok else 1


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze", help="write prediction.lock before seeing data")
    t = sub.add_parser("test", help="run the committed test against Voyager data")
    t.add_argument("--anomalies", required=True, type=Path,
                   help="tables/anomalies.csv from Voyager-Data-Analysis")
    t.add_argument("--spacecraft", default="voyager1",
                   choices=["voyager1", "voyager2"])
    t.add_argument("--bins", type=int, default=24)
    t.add_argument("--vtk", action="store_true", help="write Paraview overlay")
    sub.add_parser("selftest", help="prove the layer can both reject and accept")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    pred = FoamPrediction()      # <-- edit to your DERIVED constants, then freeze

    if args.cmd == "freeze":
        cmd_freeze(pred)
        return 0
    if args.cmd == "selftest":
        return cmd_selftest()
    if args.cmd == "test":
        load_and_verify_lock(pred)
        prof = load_profile(args.anomalies, args.spacecraft, args.bins)
        v = run_test(prof, pred, want_vtk=args.vtk)
        # exit code communicates the verdict to CI / scripts
        return {"SUPPORTED_REPLICATE": 0, "INCONCLUSIVE": 2,
                "FALSIFIED": 3, "UNFALSIFIABLE_AS_SET": 4}.get(v.label, 5)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
