"""Time-dependent snow and sea-ice conduction model (experimental).

Builds on a one-layer FTCS soil diffusion practical by Ian Renfrew.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

from ..radiation import EPS_SFC, SIGMA
from ..thermo import T_FREEZE_SW, TK0, dqsat_ice_dt, qsat_ice
from .properties import (
    BETA_K,
    C_ICE0,
    C_ICE0_V2,
    C_SNOW,
    K_ICE_FRESH,
    K_SNOW,
    L0_MU,
    RHO_ICE,
    RHO_SNOW,
    S_ICE_DEFAULT,
)
from .seb import effective_conductance, surface_albedo


class IceColumn:
    """One-dimensional snow and sea-ice heat conduction with an SEB upper boundary.

    - Fully implicit (backward Euler) time stepping on a snow layer over ice; interface
      conductances from the harmonic mean.
    - Salinity-dependent ice conductivity and heat capacity (Bitz & Lipscomb 1999, who
      cite Untersteiner 1964 for the conductivity), updated by Picard iteration each step.
    - Upper boundary: a zero-heat-capacity surface node closed by the SEB

          0 = F_cond + eps LWd + (1 - alb) SWd - eps sigma Ts^4 - SH(Ts) - LH(Ts),

      linearised about the current Ts (Newton), with SH = h_s (Ts - theta_a) and
      LH = h_l (qsat_ice(Ts) - q_a).
    - Ts capped at 0 C; the excess energy is reported as a melt flux.
    - Ice base held at ``t_base``. Initial state: the steady profile of the column's own
      layers and properties (constant conductive flux through all layers), found by
      iteration because the ice conductivity depends on temperature.

    Thickness is fixed per column: one column is one slab. :meth:`run` re-initialises the
    column after data gaps and open-water spells.

    ``v2_compat=True`` restores three v2 behaviours for the regression tests only: the
    initial profile from fresh-ice conductivity, a melt flux that could carry over from an
    earlier Picard pass, and the v2 fresh-ice heat capacity.
    """

    def __init__(
        self,
        h_ice=1.0,
        h_snow=0.2,
        n_ice=10,
        n_snow=3,
        s_ice=S_ICE_DEFAULT,
        t_base=T_FREEZE_SW,
        k_snow=K_SNOW,
        rho_snow=RHO_SNOW,
        salinity_dependent=True,
        v2_compat=False,
    ):
        self.v2_compat = bool(v2_compat)
        self.h_ice, self.h_snow = float(h_ice), float(max(h_snow, 0.0))
        self.n_snow = n_snow if self.h_snow > 0.005 else 0
        self.n_ice = int(n_ice)
        self.s_ice = float(s_ice)
        self.t_base = float(t_base)
        self.k_snow, self.rho_snow = k_snow, rho_snow
        self.salinity_dependent = salinity_dependent

        dz_s = [self.h_snow / self.n_snow] * self.n_snow if self.n_snow else []
        dz_i = [self.h_ice / self.n_ice] * self.n_ice
        self.dz = np.array(dz_s + dz_i)  # layer thickness, top down
        self.is_snow = np.array([True] * self.n_snow + [False] * self.n_ice)
        self.z = np.cumsum(self.dz) - 0.5 * self.dz  # layer-centre depth from surface
        self.J = self.dz.size
        self.T = None  # layer temperatures (C)
        self.Ts = np.nan

    # -- material properties ------------------------------------------------------------
    def _k(self, T):
        k = np.full(self.J, self.k_snow)
        if self.salinity_dependent:
            Ti = np.minimum(T[~self.is_snow], -0.1)
            k[~self.is_snow] = K_ICE_FRESH + BETA_K * self.s_ice / Ti
            k[~self.is_snow] = np.maximum(k[~self.is_snow], 0.5)
        else:
            k[~self.is_snow] = K_ICE_FRESH
        return k

    def _rhoc(self, T):
        rc = np.full(self.J, self.rho_snow * C_SNOW)
        c0 = C_ICE0_V2 if self.v2_compat else C_ICE0
        if self.salinity_dependent:
            Ti = np.minimum(T[~self.is_snow], -0.1)
            rc[~self.is_snow] = RHO_ICE * (c0 + L0_MU * self.s_ice / Ti**2)
        else:
            rc[~self.is_snow] = RHO_ICE * c0
        return rc

    def _conductances(self, T):
        """Conductances (W m-2 K-1): surface->layer 0, between layers, last layer->base."""
        k = self._k(T)
        g_top = k[0] / (0.5 * self.dz[0])
        r = 0.5 * self.dz[:-1] / k[:-1] + 0.5 * self.dz[1:] / k[1:]  # series resistance
        g_mid = 1.0 / r
        g_bot = k[-1] / (0.5 * self.dz[-1])
        return g_top, g_mid, g_bot

    # -- steady state -------------------------------------------------------------------
    def _steady(self, Ts):
        """Steady profile and conductive flux (W m-2, positive upward) for surface Ts."""
        Ts = min(float(Ts), 0.0)
        if self.v2_compat:  # v2: fresh-ice conductivity whatever the column's properties
            k_i = K_ICE_FRESH
            G = effective_conductance(self.h_ice, self.h_snow, k_i, self.k_snow)
            F = G * (self.t_base - Ts)
            T = np.empty(self.J)
            for j in range(self.J):
                zj = self.z[j]
                if self.is_snow[j]:
                    T[j] = Ts + F * zj / self.k_snow
                else:
                    T[j] = Ts + F * (self.h_snow / self.k_snow + (zj - self.h_snow) / k_i)
            return T, F
        # Constant flux F through resistances in series: surface -> layer 0 -> ... -> base.
        # The ice conductivity depends on temperature, so iterate to a fixed point.
        T = Ts + (self.t_base - Ts) * self.z / (self.h_snow + self.h_ice)
        F = 0.0
        for _ in range(100):
            g_top, g_mid, g_bot = self._conductances(T)
            r = np.concatenate(([1.0 / g_top], 1.0 / g_mid, [1.0 / g_bot]))
            F = (self.t_base - Ts) / r.sum()
            T_new = Ts + F * np.cumsum(r[:-1])
            done = np.max(np.abs(T_new - T)) < 1e-12
            T = T_new
            if done:
                break
        return T, F

    def steady_profile(self, Ts):
        """Steady layer temperatures (C) between surface Ts and the base at t_base."""
        return self._steady(Ts)[0]

    def steady_flux(self, Ts):
        """Steady conductive flux (W m-2, positive upward) for surface temperature Ts."""
        return self._steady(Ts)[1]

    def initialise(self, Ts):
        self.Ts = min(float(Ts), 0.0)
        self.T = self.steady_profile(self.Ts)

    # -- one implicit step --------------------------------------------------------------
    def step(
        self,
        dt,
        theta_a,
        q_a,
        p,
        lw_down,
        sw_down,
        h_s,
        h_l,
        albedo=None,
        n_newton=4,
        n_picard=2,
    ):
        """Advance the column by dt (s). Returns (Ts, F_cond_surface, melt_flux).

        F_cond_surface is the conductive flux reaching the surface from below
        (W m-2, positive upward).
        """
        alb = surface_albedo(self.h_snow) if albedo is None else albedo
        Tn = self.T.copy()
        Ts = self.Ts if np.isfinite(self.Ts) else -5.0
        T_new, Ts_new = Tn.copy(), Ts
        melt = 0.0
        for _ in range(n_picard):
            if not self.v2_compat:
                melt = 0.0  # a melt flux from an earlier pass must not carry over
            rhoc = self._rhoc(0.5 * (Tn + T_new))
            g_top, g_mid, g_bot = self._conductances(0.5 * (Tn + T_new))
            cap = rhoc * self.dz / dt  # W m-2 K-1
            for _ in range(n_newton):
                # unknowns x = [Ts, T_0 .. T_{J-1}]; tridiagonal, banded storage (1, 1)
                n = self.J + 1
                ab = np.zeros((3, n))
                rhs = np.zeros(n)
                # surface node: linearised SEB
                TsK = Ts_new + TK0
                lw_up0 = EPS_SFC * SIGMA * TsK**4
                dlw = 4 * EPS_SFC * SIGMA * TsK**3
                qs0 = qsat_ice(Ts_new, p)
                dqs = dqsat_ice_dt(Ts_new, p)
                # 0 = g_top (T0 - Ts) + A - lw_up0 - dlw (Ts - Ts*) - h_s (Ts - theta)
                #     - h_l (qs0 + dqs (Ts - Ts*) - q_a)
                A = EPS_SFC * lw_down + (1 - alb) * sw_down
                diag_s = -(g_top + dlw + h_s + h_l * dqs)
                rhs_s = -(
                    A - lw_up0 + dlw * Ts_new + h_s * theta_a - h_l * (qs0 - dqs * Ts_new - q_a)
                )
                ab[1, 0] = diag_s
                ab[0, 1] = g_top  # coefficient of T0 in row 0
                rhs[0] = rhs_s
                # layer rows
                for j in range(self.J):
                    row = j + 1
                    g_up = g_top if j == 0 else g_mid[j - 1]
                    g_dn = g_bot if j == self.J - 1 else g_mid[j]
                    ab[1, row] = -(cap[j] + g_up + g_dn)
                    ab[2, row - 1] = g_up  # coefficient of node above
                    if j < self.J - 1:
                        ab[0, row + 1] = g_dn  # coefficient of node below
                    rhs[row] = -cap[j] * Tn[j] - (g_bot * self.t_base if j == self.J - 1 else 0.0)
                x = solve_banded((1, 1), ab, rhs)
                Ts_new, T_new = x[0], x[1:]
                if Ts_new > 0.0:  # melting surface: fix Ts = 0
                    Ts_new = 0.0
                    ab2, rhs2 = ab[:, 1:].copy(), rhs[1:].copy()
                    rhs2[0] -= g_top * 0.0  # Ts = 0 contributes g_top * 0
                    ab2[2, -1] = 0.0
                    T_new = solve_banded((1, 1), ab2, rhs2)
                    F_c = g_top * (T_new[0] - 0.0)
                    melt = (
                        A
                        - EPS_SFC * SIGMA * TK0**4
                        - h_s * (0.0 - theta_a)
                        - h_l * (qsat_ice(0.0, p) - q_a)
                        + F_c
                    )
                    break
        self.T, self.Ts = T_new, Ts_new
        F_c = g_top * (self.T[0] - self.Ts)
        return self.Ts, F_c, max(melt, 0.0)

    # -- drive along a time series ------------------------------------------------------
    def run(
        self,
        time_s,
        theta_a,
        q_a,
        p,
        lw_down,
        sw_down,
        h_s,
        h_l,
        Ts_init,
        sic=None,
        sic_reset=0.15,
        max_gap_s=3 * 3600.0,
        albedo=None,
    ):
        """Integrate along a time series (all arrays length N; time_s in seconds).

        The column is (re)initialised to the steady profile for Ts_init[i] at the first
        valid step, after a gap longer than max_gap_s, and on re-entering ice after
        SIC < sic_reset. Steps with any non-finite input, or out of ice, are skipped
        (outputs NaN).

        Returns a dict of arrays: Ts, F_cond, melt, reinit (bool).
        """
        N = len(time_s)
        out = {k: np.full(N, np.nan) for k in ("Ts", "F_cond", "melt")}
        out["reinit"] = np.zeros(N, dtype=bool)
        last_t = None
        inputs = [theta_a, q_a, p, lw_down, sw_down, h_s, h_l, Ts_init]
        for i in range(N):
            vals = [np.asarray(a)[i] for a in inputs]
            ok = all(np.isfinite(v) for v in vals)
            in_ice = True if sic is None else (np.isfinite(sic[i]) and sic[i] >= sic_reset)
            if not (ok and in_ice):
                last_t = None if not in_ice else last_t
                continue
            if last_t is None or (time_s[i] - last_t) > max_gap_s or self.T is None:
                self.initialise(vals[7])
                out["reinit"][i] = True
                out["Ts"][i] = self.Ts
                if self.v2_compat:  # v2: default fresh-ice and default-snow conductance
                    out["F_cond"][i] = effective_conductance(self.h_ice, self.h_snow) * (
                        self.t_base - self.Ts
                    )
                else:
                    out["F_cond"][i] = self.steady_flux(self.Ts)
                out["melt"][i] = 0.0
            else:
                Ts, Fc, m = self.step(time_s[i] - last_t, *vals[:7], albedo=albedo)
                out["Ts"][i], out["F_cond"][i], out["melt"][i] = Ts, Fc, m
            last_t = time_s[i]
        return out
