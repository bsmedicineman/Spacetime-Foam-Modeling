#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
solar_foam_animation.py

Create a time‑dependent foam‑distortion field produced by the Sun,
planets and a handful of large moons.  The baseline foam grid is
imported from the program that you already have (baseline_foam_map.py).

Output format: VTK Structured Grid (vts) + a tiny JSON side‑car that
contains the body positions, radii and display colour for the current
frame.  Any modern open‑source visualiser (Paraview, Blender, VisIt,
vtk.js) can read this pair and animate a “living” solar system where
foam waves are coloured according to amplitude.

All physics follows the equations in **Master Equation Reference.docx**:

* 9.1 – Gravitational potential from the foam field  (ϕ_g = α_g * M / r)
* 9.3 – Modified Poisson (source term = Σ M_i / r_i) → we use the
  linear superposition of the potentials.
* 3.4 – Warp‑velocity from strain  (v_w = ε / L)  → we turn the potential
  gradient into a strain field ε = ∇ϕ_g / k_s, then into a velocity‑like
  “foam‑wave” amplitude.
* 4.3 – Force‑balance term (radiation pressure) – we add a tiny EM
  contribution for the Sun (P ≈ 0.173 N, B ≈ 0.3 N, S ≈ 0.1 N from the
  reference).  The EM term is treated as an extra scalar source
  ϕ_EM = β_EM * P / r².
* 5.17 – Foam‑field energy (quadratic) – we compute a simple proxy
  `phi_dist = sqrt(ϕ_g² + ϕ_EM²)` that the visualiser colours.

The code is deliberately **explicit** (no hidden magic) so you can
swap any block for a more elaborate model later.

Author:  BSMEDICINEMAN (6.13.2026)   License: MIT
"""

# ----------------------------------------------------------------------
# 0️⃣  IMPORTS & GLOBAL SETTINGS
# ----------------------------------------------------------------------
import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm
from pyvtk import StructuredGrid, PointData, VtkData, Scalars, Vectors

# ----------------------------------------------------------------------
# 1️⃣  PHYSICAL CONSTANTS (taken from the reference)
# ----------------------------------------------------------------------
AU_M = 1.495978707e11          # metres
G   = 6.67430e-11              # m³ kg⁻¹ s⁻²
C   = 299_792_458.0            # m s⁻¹

# Foam‑related constants (see §3.2, §9.1)
ALPHA_G = 1.0e-27   # coupling constant for foam‑gravity (ϕ_g = α_g * M / r)   – placeholder
K_STR   = 1.0e-6    # strain‑to‑gradient conversion factor (ε = |∇ϕ_g| / K_STR) – placeholder
BETA_EM = 1.0e-12   # EM‑foam coupling (ϕ_EM = β_EM * P / r²)                – placeholder

# ----------------------------------------------------------------------
# 2️⃣  DATA STRUCTURES
# ----------------------------------------------------------------------
@dataclass
class Body:
    """Simple container for a gravitating body."""
    name: str
    mass_kg: float
    radius_m: float          # physical radius (used for visual size)
    color_rgb: Tuple[int, int, int]   # 0‑255 for visualiser
    # Keplerian elements (all angles in radians, epoch J2000)
    a_au: float              # semi‑major axis
    e: float                 # eccentricity
    i_deg: float             # inclination
    Ω_deg: float             # longitude of ascending node
    ω_deg: float             # argument of periapsis
    M0_deg: float            # mean anomaly at epoch

    # derived quantities (filled after init)
    mu: float = 0.0          # GM (m³ s⁻²)
    a_m: float = 0.0
    i: float = 0.0
    Ω: float = 0.0
    ω: float = 0.0
    M0: float = 0.0

    def __post_init__(self):
        self.mu = G * self.mass_kg
        self.a_m = self.a_au * AU_M
        self.i = math.radians(self.i_deg)
        self.Ω = math.radians(self.Ω_deg)
        self.ω = math.radians(self.ω_deg)
        self.M0 = math.radians(self.M0_deg)


# ----------------------------------------------------------------------
# 3️⃣  EPHEMERIS – SIMPLE KEPLERIAN PROPAGATOR
# ----------------------------------------------------------------------
def keplerian_position(body: Body, t_seconds: float) -> np.ndarray:
    """
    Return the heliocentric Cartesian position (m) of *body* at elapsed
    time `t_seconds` after the J2000 epoch, using the classic
    Keplerian solution (mean motion n = sqrt(μ/a³)).

    This is sufficient for a visual sandbox.  For high‑precision work
    replace this function with a SPICE call (e.g. spiceypy) or a
    numerical integrator.
    """
    # 1️⃣ mean motion
    n = math.sqrt(body.mu / body.a_m ** 3)                     # rad/s

    # 2️⃣ mean anomaly at time t
    M = body.M0 + n * t_seconds                                 # rad
    M = M % (2 * math.pi)

    # 3️⃣ solve Kepler's equation (Newton‑Raphson, 5 iterations)
    E = M
    for _ in range(5):
        E = E - (E - body.e * math.sin(E) - M) / (1 - body.e * math.cos(E))

    # 4️⃣ true anomaly
    ν = 2 * math.atan2(math.sqrt(1 + body.e) * math.sin(E / 2),
                      math.sqrt(1 - body.e) * math.cos(E / 2))

    # 5️⃣ distance from focus
    r = body.a_m * (1 - body.e * math.cos(E))

    # 6️⃣ orbital plane coordinates
    x_orb = r * math.cos(ν)
    y_orb = r * math.sin(ν)
    z_orb = 0.0

    # 7️⃣ rotate to ecliptic J2000 frame
    # Rotation matrix: Rz(-Ω) * Rx(-i) * Rz(-ω)
    cosΩ, sinΩ = math.cos(-body.Ω), math.sin(-body.Ω)
    cosi, sini = math.cos(-body.i), math.sin(-body.i)
    cosω, sinω = math.cos(-body.ω), math.sin(-body.ω)

    # Build the 3×3 matrix (hard‑coded for speed)
    R = np.array([
        [cosΩ*cosω - sinΩ*sinω*cosi,
         -cosΩ*sinω - sinΩ*cosω*cosi,
         sinΩ*sini],
        [sinΩ*cosω + cosΩ*sinω*cosi,
         -sinΩ*sinω + cosΩ*cosω*cosi,
         -cosΩ*sini],
        [sinω*sini,
         cosω*sini,
         cosi]
    ])

    pos = R @ np.array([x_orb, y_orb, z_orb])
    return pos  # metres, heliocentric


# ----------------------------------------------------------------------
# 4️⃣  LOAD THE BASELINE FOAM GRID
# ----------------------------------------------------------------------
def load_baseline_grid() -> Tuple[np.ndarray, float, float]:
    """
    Import the FoamMap class from the baseline program, instantiate it,
    and return the three‑dimensional coordinate arrays together with the
    uniform cell spacing (in metres).

    The returned tuple is:
        (grid_xyz, spacing_m, origin_m)

    *grid_xyz* is a 3‑tuple of (X, Y, Z) arrays shaped (nx, ny, nz) that
    contain the physical coordinates of each cell centre.
    """
    # The baseline program lives in the same directory – we import it
    import baseline_foam_map as bfm
    foam = bfm.FoamMap()
    # The coordinate arrays are already public members
    return (foam.x_grid, foam.y_grid, foam.z_grid), foam.spacing_m, -foam.r_max_m


# ----------------------------------------------------------------------
# 5️⃣  PHYSICS – FOAM DISTORTION FROM A SINGLE BODY
# ----------------------------------------------------------------------
def body_foam_contribution(body: Body,
                          X: np.ndarray,
                          Y: np.ndarray,
                          Z: np.ndarray,
                          pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the scalar foam disturbance (ϕ_dist) and the foam‑gradient
    force vector (F_foam = -∇ϕ_dist) produced by *body* at the current
    position `pos` (3‑vector, metres).

    The model follows the equations in the reference:

    * 9.1 – ϕ_g = α_g * M / r
    * 9.3 – we simply superpose contributions from all bodies.
    * 3.4 – strain ε = |∇ϕ_g| / K_STR
    * 4.3 – add a tiny electromagnetic term ϕ_EM = β_EM * P / r²
      (P is the radiation‑pressure force listed in §4.3; we use the
      Sun’s baseline value P≈0.173 N).

    Returned arrays have the same shape as the input grid.
    """
    # Vector from body centre to every grid point
    dx = X - pos[0]
    dy = Y - pos[1]
    dz = Z - pos[2]

    # Distance (avoid divide‑by‑zero at the body centre)
    r = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-12

    # ---- 1️⃣ Gravitational foam scalar ----------
    phi_g = ALPHA_G * body.mass_kg / r                     # dimensionless

    # ---- 2️⃣ EM foam (only for the Sun, others are negligible) ----------
    if body.name.lower() == "sun":
        # Use the reference radiation‑pressure value P≈0.173 N (Sec. 4.3)
        P_sun = 0.173
        phi_em = BETA_EM * P_sun / (r * r)                # falls as 1/r²
    else:
        phi_em = 0.0

    # ---- 3️⃣ Total disturbance (scalar) ----------
    phi_dist = np.sqrt(phi_g * phi_g + phi_em * phi_em)   # Eq. 5.17 proxy

    # ---- 4️⃣ Foam‑gradient force (vector) ----------
    # ∇ϕ_g = -α_g * M * r̂ / r²   ; we reuse phi_g to avoid recompute
    grad_factor = -ALPHA_G * body.mass_kg / (r * r * r)   # = -α_g M / r³
    Fx = grad_factor * dx
    Fy = grad_factor * dy
    Fz = grad_factor * dz

    # Turn the gradient into a *force*‑like vector (the sign is already
    # correct for a restoring attraction toward the mass).
    foam_force = np.stack([Fx, Fy, Fz], axis=-1)          # shape (nx,ny,nz,3)

    return phi_dist, foam_force


# ----------------------------------------------------------------------
# 6️⃣  COMBINE ALL BODIES FOR ONE TIME‑STEP
# ----------------------------------------------------------------------
def compute_foam_step(bodies: List[Body],
                      X: np.ndarray,
                      Y: np.ndarray,
                      Z: np.ndarray,
                      t_seconds: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return the total scalar distortion field `phi_total` and the
    vector field `foam_force_total` for the whole grid at simulation
    time `t_seconds`.  The algorithm:

    1. Loop over every body,
    2. Get its heliocentric position,
    3. Add its contribution (scalar + vector) to the accumulator.
    """
    # Initialise accumulators with zeros of the proper shape
    shape = X.shape
    phi_total = np.zeros(shape, dtype=np.float64)
    force_total = np.zeros(shape + (3,), dtype=np.float64)

    for body in bodies:
        pos = keplerian_position(body, t_seconds)          # metres
        phi_body, force_body = body_foam_contribution(body, X, Y, Z, pos)

        phi_total   += phi_body
        force_total += force_body

    return phi_total, force_total


# ----------------------------------------------------------------------
# 7️⃣  WRITE VTK + JSON SIDE‑CAR FOR ONE FRAME
# ----------------------------------------------------------------------
def write_frame(frame_idx: int,
                phi: np.ndarray,
                force: np.ndarray,
                bodies: List[Body],
                grid_origin: float,
                spacing: float,
                out_dir: Path):
    """
    *phi*   – scalar foam distortion (dimensionless)  
    *force* – vector foam‑gradient force (same units as φ‑gradient, i.e.
              dimensionless m⁻¹)

    The VTK file contains two data arrays:
        * "phi_dist"  (scalar) – colour‑mapped by the visualiser
        * "foam_force" (vector) – optional glyph/arrow display

    The JSON side‑car stores a list of bodies with
        - name
        - position (AU)
        - radius (km, for scaling in the visualiser)
        - colour (RGB 0‑255)
    """
    nx, ny, nz = phi.shape
    # Build the StructuredGrid object (origin is the lower‑left‑back corner)
    grid = StructuredGrid(
        (nx, ny, nz),
        spacing=(spacing, spacing, spacing),
        origin=(grid_origin, grid_origin, grid_origin)
    )

    # Attach point data
    pdata = PointData(
        Scalars(phi, name="phi_dist"),
        Vectors(force, name="foam_force")
    )
    vtk = VtkData(grid, pdata)

    vtk_path = out_dir / f"foam_{frame_idx:05d}.vts"
    vtk.tofile(str(vtk_path), format="binary")

    # ---- JSON side‑car ----
    bodies_json = []
    for body in bodies:
        pos = keplerian_position(body, frame_idx * dt_seconds)  # metres
        bodies_json.append({
            "name": body.name,
            "pos_au": [pos[0] / AU_M, pos[1] / AU_M, pos[2] / AU_M],
            "radius_km": body.radius_m / 1e3,
            "color_rgb": body.color_rgb
        })

    json_path = out_dir / f"bodies_{frame_idx:05d}.json"
    with open(json_path, "w") as f:
        json.dump(bodies_json, f, indent=2)

# ----------------------------------------------------------------------
# 8️⃣  BUILD THE SOLAR‑SYSTEM BODY LIST
# ----------------------------------------------------------------------
def build_solar_system_bodies() -> List[Body]:
    """
    Returns a list of Body objects for the Sun, 8 planets,
    Pluto (dwarf), and a few large moons (Moon, Europa, Ganymede,
    Titan).  The orbital elements are taken from J2000 NASA data
    (public domain).  Colours are approximate visual colours.
    """
    # Helper to convert km → m for radii
    km = 1e3

    bodies = [
        Body(
            name="Sun",
            mass_kg=1.9885e30,
            radius_m=696_340 * km,
            color_rgb=(255, 204, 0),
            a_au=0.0,
            e=0.0,
            i_deg=0.0,
            Ω_deg=0.0,
            ω_deg=0.0,
            M0_deg=0.0
        ),
        Body(
            name="Mercury",
            mass_kg=3.3011e23,
            radius_m=2_439.7 * km,
            color_rgb=(180, 180, 180),
            a_au=0.387098,
            e=0.205630,
            i_deg=7.00487,
            Ω_deg=48.33167,
            ω_deg=29.12478,
            M0_deg=174.796
        ),
        Body(
            name="Venus",
            mass_kg=4.8675e24,
            radius_m=6_051.8 * km,
            color_rgb=(255, 200, 150),
            a_au=0.723332,
            e=0.006772,
            i_deg=3.39471,
            Ω_deg=76.68069,
            ω_deg=54.85229,
            M0_deg=50.115
        ),
        Body(
            name="Earth",
            mass_kg=5.97237e24,
            radius_m=6_371.0 * km,
            color_rgb=(0, 102, 255),
            a_au=1.000000,
            e=0.0167086,
            i_deg=0.00005,
            Ω_deg=-11.26064,
            ω_deg=102.93735,
            M0_deg=357.51716
        ),
        Body(
            name="Mars",
            mass_kg=6.4171e23,
            radius_m=3_389.5 * km,
            color_rgb=(210, 105, 30),
            a_au=1.523679,
            e=0.0934,
            i_deg=1.850,
            Ω_deg=49.558,
            ω_deg=286.502,
            M0_deg=19.3564
        ),
        Body(
            name="Jupiter",
            mass_kg=1.8982e27,
            radius_m=69_911 * km,
            color_rgb=(210, 180, 140),
            a_au=5.204267,
            e=0.0489,
            i_deg=1.303,
            Ω_deg=100.464,
            ω_deg=273.867,
            M0_deg=20.020
        ),
        Body(
            name="Saturn",
            mass_kg=5.6834e26,
            radius_m=58_232 * km,
            color_rgb=(210, 190, 140),
            a_au=9.582017,
            e=0.0565,
            i_deg=2.485,
            Ω_deg=113.665,
            ω_deg=339.392,
            M0_deg=317.020
        ),
        Body(
            name="Uranus",
            mass_kg=8.6810e25,
            radius_m=25_362 * km,
            color_rgb=(150, 220, 255),
            a_au=19.189165,
            e=0.0472,
            i_deg=0.773,
            Ω_deg=74.006,
            ω_deg=96.998857,
            M0_deg=142.238
        ),
        Body(
            name="Neptune",
            mass_kg=1.02413e26,
            radius_m=24_622 * km,
            color_rgb=(100, 150, 255),
            a_au=30.069922,
            e=0.0086,
            i_deg=1.770,
            Ω_deg=131.784,
            ω_deg=272.846,
            M0_deg=256.228
        ),
        Body(
            name="Pluto",
            mass_kg=1.303e22,
            radius_m=1_188.3 * km,
            color_rgb=(200, 200, 255),
            a_au=39.482,
            e=0.2488,
            i_deg=17.16,
            Ω_deg=110.30,
            ω_deg=113.78,
            M0_deg=14.53
        ),
        # ----- Large moons (treated as separate bodies, orbiting their planet) -----
        Body(
            name="Moon",
            mass_kg=7.342e22,
            radius_m=1_737.1 * km,
            color_rgb=(200, 200, 200),
            a_au=0.00257,          # 384 400 km ≈ 0.00257 AU from Earth
            e=0.0549,
            i_deg=5.145,
            Ω_deg=0.0,
            ω_deg=0.0,
            M0_deg=0.0
        ),
        Body(
            name="Europa",
            mass_kg=4.799e22,
            radius_m=1_560.8 * km,
            color_rgb=(180, 220, 255),
            a_au=0.000083,        # 670 900 km ≈ 0.000083 AU from Jupiter
            e=0.009,
            i_deg=0.470,
            Ω_deg=0.0,
            ω_deg=0.0,
            M0_deg=0.0
        ),
        Body(
            name="Ganymede",
            mass_kg=1.4819e23,
            radius_m=2_634.1 * km,
            color_rgb=(150, 190, 255),
            a_au=0.000152,        # 1 070 400 km ≈ 0.000152 AU
            e=0.0013,
            i_deg=0.177,
            Ω_deg=0.0,
            ω_deg=0.0,
            M0_deg=0.0
        ),
        Body(
            name="Titan",
            mass_kg=1.3452e23,
            radius_m=2_575.5 * km,
            color_rgb=(255, 210, 180),
            a_au=0.000886,        # 1 221 830 km ≈ 0.000886 AU from Saturn
            e=0.0288,
            i_deg=0.348,
            Ω_deg=0.0,
            ω_deg=0.0,
            M0_deg=0.0
        ),
    ]

    return bodies


# ----------------------------------------------------------------------
# 9️⃣  MAIN DRIVER
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a time‑dependent foam‑distortion field "
                    "for the Solar System and export it as VTK + JSON.")
    parser.add_argument("--duration-days", type=float, default=365,
                        help="Total simulated time in days (default: 365).")
    parser.add_argument("--dt-hours", type=float, default=1,
                        help="Time step in hours (default: 1 h).")
    parser.add_argument("--output", type=Path, default=Path("./foam_frames"),
                        help="Directory where VTK/JSON frames will be written.")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 9.1  Prepare output folder
    # ------------------------------------------------------------------
    args.output.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 9.2  Load baseline grid
    # ------------------------------------------------------------------
    (X, Y, Z), spacing_m, origin_m = load_baseline_grid()
    nx, ny, nz = X.shape
    print(f"[Info] Baseline grid loaded: {nx}×{ny}×{nz} cells, "
          f"spacing = {spacing_m:.2e} m")

    # ------------------------------------------------------------------
    # 9.3  Build body list
    # ------------------------------------------------------------------
    bodies = build_solar_system_bodies()
    print(f"[Info] {len(bodies)} bodies loaded (including moons).")

    # ------------------------------------------------------------------
    # 9.4  Time stepping
    # ------------------------------------------------------------------
    total_seconds = args.duration_days * 24 * 3600
    dt_seconds = args.dt_hours * 3600
    n_steps = int(np.ceil(total_seconds / dt_seconds))

    print(f"[Info] Running {n_steps} steps (Δt = {args.dt_hours:.2f} h).")

    for step in tqdm(range(n_steps), desc="Generating frames"):
        t = step * dt_seconds

        # ---- Compute foam field for this instant ----
        phi_tot, force_tot = compute_foam_step(bodies, X, Y, Z, t)

        # ---- Write VTK + JSON side‑car ----
        write_frame(step,
                    phi_tot,
                    force_tot,
                    bodies,
                    origin_m,
                    spacing_m,
                    args.output)

    print(f"\n[Done] {n_steps} VTK frames written to {args.output.resolve()}")
    print("Load the *.vts files in Paraview (or any VTK‑compatible viewer).")
    print("The accompanying JSON files contain the planet positions and colours.")
