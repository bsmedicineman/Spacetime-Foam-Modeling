#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
baseline_foam_map.py

Create a baseline spacetime‑foam matrix for the Solar System out to
10 AU beyond the Sun (≈ 20 AU per side).  The map is a sandbox that
future programs can query for craft positioning and can stream to a
3‑D visualiser.

All equations are taken from **Master Equation Reference.docx**:

* §3.1 – scalar foam field ϕ (baseline = 0)
* §3.2 – foam fluid properties (ρ_f, σ_f, η_f)
* §6.1 – temporal refractive index n_t  (baseline = 1)
* §3.4 – warp‑velocity from strain (unused here because strain = 0)
* §4.3 – force‑balance (baseline gives zero net force)
* §5.17 – foam field energy (quadratic ≈ 0 for the undisturbed case)

The map assumes **no disturbances**; therefore strain, curvature,
and any derived “distortion” fields are identically zero.
Future models can replace the placeholder values with measured
perturbations.

Author:  BSMEDICINEMAN (06.13.2026)   License: MIT
"""

import json
import socket
import struct
import threading
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

# ----------------------------------------------------------------------
# 1️⃣  GLOBAL CONSTANTS  (taken from the reference)
# ----------------------------------------------------------------------
AU_M = 1.495978707e11                     # 1 AU in metres
R_MAX_AU = 20.0                           # cube side = 20 AU (±10 AU from Sun)
GRID_SPACING_AU = 0.1                     # resolution – 0.1 AU ≈ 15 000 km
GRID_SPACING_M = GRID_SPACING_AU * AU_M

# Baseline foam fluid properties – see §3.2
# (the reference gives symbolic placeholders; we use plausible order‑of‑magnitude values)
RHO_F = 1e-21          # kg m⁻³   – effective spacetime density
SIGMA_F = 1e-5         # N m⁻¹   – surface tension of the foam
ETA_F = 1e-30          # Pa·s    – dynamic viscosity

# Baseline scalar foam field ϕ – see §3.1 (small‑perturbation form)
PHI_0 = 0.0            # dimensionless; undisturbed foam = 0

# Temporal refractive index – see §6.1 (n_t = 1 for flat spacetime)
N_T = 1.0

# Derived constants (used for optional diagnostics)
#   Reciprocal energy–time identification (§1.10)
#   In natural units c = ħ = 1, so we keep them as 1 for the baseline.
C_LIGHT = 299_792_458.0    # m s⁻¹
H_BAR   = 1.054_571_8e-34 # J s

# ----------------------------------------------------------------------
# 2️⃣  GRID CONSTRUCTION
# ----------------------------------------------------------------------
class FoamMap:
    """
    Holds a 3‑D numpy array with the baseline foam quantities.
    The grid spans [-10 AU, +10 AU] in each Cartesian direction.
    """

    def __init__(self,
                 r_max_au: float = R_MAX_AU / 2,
                 spacing_au: float = GRID_SPACING_AU):
        self.r_max_m = r_max_au * AU_M
        self.spacing_m = spacing_au * AU_M

        # Number of points per axis (including both ends)
        self.nx = self.ny = self.nz = int(np.round(2 * self.r_max_m / self.spacing_m)) + 1

        # Create coordinate arrays (centers of cells)
        axis = np.linspace(-self.r_max_m,
                           +self.r_max_m,
                           self.nx,
                           endpoint=True)
        self.x_grid, self.y_grid, self.z_grid = np.meshgrid(axis,
                                                            axis,
                                                            axis,
                                                            indexing='ij')

        # Allocate the baseline fields (all uniform, but stored for easy slicing)
        self.density = np.full((self.nx, self.ny, self.nz), RHO_F, dtype=np.float64)
        self.surface_tension = np.full_like(self.density, SIGMA_F)
        self.viscosity = np.full_like(self.density, ETA_F)
        self.scalar_foam = np.full_like(self.density, PHI_0)          # ϕ
        self.temporal_index = np.full_like(self.density, N_T)        # n_t

        # For completeness we also store a “distortion” field that is zero now
        self.strain = np.zeros_like(self.density)                    # ε (from §3.4)
        self.foam_energy = np.zeros_like(self.density)               # ℰ (from §5.17)

        # Pre‑compute a lookup dictionary for fast per‑coordinate queries
        self._build_lookup()

    # ------------------------------------------------------------------
    # 2️⃣.1  Helper: convert physical (m) → grid index
    # ------------------------------------------------------------------
    def _pos_to_idx(self, x: float, y: float, z: float) -> Tuple[int, int, int]:
        """Return integer grid indices for a point given in metres."""
        ix = int(round((x + self.r_max_m) / self.spacing_m))
        iy = int(round((y + self.r_max_m) / self.spacing_m))
        iz = int(round((z + self.r_max_m) / self.spacing_m))

        # Clamp to array bounds (out‑of‑range points are considered “outside” the map)
        ix = np.clip(ix, 0, self.nx - 1)
        iy = np.clip(iy, 0, self.ny - 1)
        iz = np.clip(iz, 0, self.nz - 1)
        return ix, iy, iz

    # ------------------------------------------------------------------
    # 2️⃣.2  Build a tiny cache for fast look‑ups (optional)
    # ------------------------------------------------------------------
    def _build_lookup(self):
        """Create a flat view that can be indexed with a 1‑D integer."""
        self._flat_density = self.density.ravel()
        self._flat_surface = self.surface_tension.ravel()
        self._flat_visc = self.viscosity.ravel()
        self._flat_phi = self.scalar_foam.ravel()
        self._flat_nt = self.temporal_index.ravel()
        self._flat_strain = self.strain.ravel()
        self._flat_energy = self.foam_energy.ravel()
        self._shape = self.density.shape
        self._stride_y = self._shape[2]
        self._stride_x = self._shape[1] * self._shape[2]

    # ------------------------------------------------------------------
    # 2️⃣.3  Public query API
    # ------------------------------------------------------------------
    def get_foam_at(self, x: float, y: float, z: float) -> Dict[str, float]:
        """
        Return a dictionary with the baseline foam quantities at the
        supplied Cartesian coordinates (in metres).

        The returned values are *exactly* the baseline (no distortion)
        because this map represents the undisturbed state.
        """
        ix, iy, iz = self._pos_to_idx(x, y, z)
        flat_index = ix * self._stride_x + iy * self._stride_y + iz

        return {
            "rho_f":   float(self._flat_density[flat_index]),   # kg/m³
            "sigma_f": float(self._flat_surface[flat_index]),   # N/m
            "eta_f":   float(self._flat_visc[flat_index]),     # Pa·s
            "phi":    float(self._flat_phi[flat_index]),       # scalar foam (dimensionless)
            "n_t":    float(self._flat_nt[flat_index]),        # temporal refractive index
            "strain": float(self._flat_strain[flat_index]),    # ε (zero)
            "energy": float(self._flat_energy[flat_index]),    # ℰ (zero)
        }

    # ------------------------------------------------------------------
    # 2️⃣.4  Slice‑by‑slice iterator (used by the streaming server)
    # ------------------------------------------------------------------
    def iter_slices(self):
        """
        Yield each XY‑plane (constant Z) as a JSON‑serialisable dict.
        This keeps the memory footprint modest for the receiver.
        """
        for iz in range(self.nz):
            slice_dict = {
                "z_index": iz,
                "z_m":   float(self.z_grid[0, 0, iz]),
                "rho_f": self.density[:, :, iz].tolist(),
                "sigma_f": self.surface_tension[:, :, iz].tolist(),
                "eta_f": self.viscosity[:, :, iz].tolist(),
                "phi": self.scalar_foam[:, :, iz].tolist(),
                "n_t": self.temporal_index[:, :, iz].tolist(),
                # strain and energy are all zero – omit to save bandwidth
            }
            yield slice_dict

# ----------------------------------------------------------------------
# 3️⃣  STREAMING SERVER
# ----------------------------------------------------------------------
class FoamStreamingServer(threading.Thread):
    """
    Tiny TCP server that streams the baseline foam map slice‑by‑slice.
    The protocol is **JSON‑lines** (one JSON object per line).  A client
    can connect, receive the header, then read each slice until the
    line `"END"` is encountered.

    Example client (pseudo‑code):
        s = socket.create_connection(("localhost", 7777))
        for line in s.makefile():
            if line.strip() == b'END': break
            slice = json.loads(line)
            render(slice)
    """

    def __init__(self, foam_map: FoamMap, host: str = "127.0.0.1", port: int = 7777):
        super().__init__(daemon=True)
        self.foam_map = foam_map
        self.host = host
        self.port = port
        self._shutdown = threading.Event()

    def run(self):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as srv:
            srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            srv.bind((self.host, self.port))
            srv.listen(1)
            print(f"[FoamServer] Listening on {self.host}:{self.port}")

            while not self._shutdown.is_set():
                try:
                    srv.settimeout(1.0)
                    conn, addr = srv.accept()
                except socket.timeout:
                    continue
                print(f"[FoamServer] Connection from {addr}")
                with conn:
                    self._handle_client(conn)

    def _handle_client(self, conn: socket.socket):
        """
        Send a simple JSON header, then stream each slice.
        """
        # Header – includes grid meta‑data useful for the visualiser
        header = {
            "type": "foam_map_header",
            "grid": {
                "nx": self.foam_map.nx,
                "ny": self.foam_map.ny,
                "nz": self.foam_map.nz,
                "spacing_m": self.foam_map.spacing_m,
                "origin_m": [-self.foam_map.r_max_m,
                             -self.foam_map.r_max_m,
                             -self.foam_map.r_max_m],
            },
            "units": {
                "density": "kg/m^3",
                "surface_tension": "N/m",
                "viscosity": "Pa·s",
                "temporal_index": "dimensionless",
            },
            "baseline": {
                "phi": PHI_0,
                "n_t": N_T,
            }
        }
        conn.sendall((json.dumps(header) + "\n").encode("utf-8"))

        # Stream each XY‑plane
        for slice_dict in self.foam_map.iter_slices():
            line = json.dumps(slice_dict) + "\n"
            conn.sendall(line.encode("utf-8"))

        # Terminator
        conn.sendall(b"END\n")
        print("[FoamServer] Finished streaming to client")

    def stop(self):
        self._shutdown.set()


# ----------------------------------------------------------------------
# 4️⃣  MAIN – build map, start server, provide a tiny demo query
# ----------------------------------------------------------------------
def main():
    # 4.1  Build the baseline foam matrix
    foam = FoamMap()
    print(f"[Main] Foam map built: {foam.nx}×{foam.ny}×{foam.nz} cells "
          f"({foam.nx * foam.ny * foam.nz:,} total)")

    # 4.2  Start the streaming server in a background thread
    server = FoamStreamingServer(foam)
    server.start()
    print("[Main] Streaming server started (Ctrl‑C to stop)")

    # 4.3  Simple interactive demo – query a few positions
    #      (positions are given in AU for readability)
    demo_points = [
        (0.0, 0.0, 0.0),          # Sun centre
        (1.0, 0.0, 0.0),          # 1 AU (Earth orbit)
        (5.2, 0.0, 0.0),          # Jupiter orbit
        (10.0, 0.0, 0.0),         # Edge of inner Solar System
        (15.0, 0.0, 0.0),         # 5 AU into interstellar space
    ]

    for au_x, au_y, au_z in demo_points:
        x, y, z = au_x * AU_M, au_y * AU_M, au_z * AU_M
        props = foam.get_foam_at(x, y, z)
        print(f"\n[Demo] Position ({au_x:.1f}, {au_y:.1f}, {au_z:.1f}) AU → "
              f"ρ_f={props['rho_f']:.2e} kg/m³, σ_f={props['sigma_f']:.2e} N/m, "
              f"η_f={props['eta_f']:.2e} Pa·s, n_t={props['n_t']:.1f}")

    # Keep the process alive until the user hits Ctrl‑C
    try:
        while True:
            pass
    except KeyboardInterrupt:
        print("\n[Main] Shutting down...")
        server.stop()
        server.join()
        print("[Main] Bye!")


if __name__ == "__main__":
    main()
