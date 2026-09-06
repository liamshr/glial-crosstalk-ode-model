"""Reusable rate equations and solve_ivp wrapper for the nine-state model."""

import numpy as np
from scipy.integrate import solve_ivp

from .parameters import DEFAULT_INITIAL_CONDITIONS, DEFAULT_PARAMETERS, ModelParameters

STATE_NAMES = ("M1", "M2", "IL1", "IL12", "IL10", "IL4", "Aq", "Ap", "D")


def _hill(x, threshold, exponent):
    x = max(float(x), 0.0)
    numerator = x**exponent
    return numerator / (max(threshold, 0.0) ** exponent + numerator)


def rates(y, params=DEFAULT_PARAMETERS):
    """Return auxiliary rates used by the ODE RHS."""
    M1, M2, IL1, IL12, IL10, IL4, Aq, Ap, D = np.maximum(y, 0.0)
    rm1 = _hill(params.k_n1 * IL1 + params.k_n12 * IL12, params.b_n, params.x_n)
    rm1 *= (1 + params.beta_Ap_M1 * Ap) / (1 + params.beta_Aq_inhib_M1 * Aq)
    rm1 /= 1 + ((IL10 + IL4) / params.a_inf1) ** 2
    rm2 = _hill(params.k_n4 * IL4 + params.k_n10 * IL10, params.y_n, params.z_n)
    rm2 *= 1 + params.beta_Aq_M2 * Aq
    rms = _hill(params.tau_n4 * IL4 + params.tau_n10 * IL10, params.m_n, params.g_n)
    mr = params.s_mr / (rm1 + rm2 + params.mu_mr)
    rp = M1 * (params.k_M1base + _hill(IL12 + params.kcd * D, params.v_n, params.h_n))
    rp /= 1 + ((IL10 + IL4) / params.a_inf1) ** 2
    rt = (params.k_tbase + _hill(IL4 + params.k_tn12 * IL12, params.r_n, params.c_n))
    rt /= 1 + (IL10 / params.a_inf2) ** 2
    ra = M2 * (params.k_M2base + _hill(params.kc4 * IL4, params.w_n, params.q_n))
    ra /= 1 + (IL10 / params.a_inf2) ** 2
    return {"Rm1": rm1, "Rm2": rm2, "Rms": rms, "mr": mr, "Rp": rp, "Rt": rt, "Ra": ra}


def rhs(t, y, params=DEFAULT_PARAMETERS):
    """Compute dY/dt in STATE_NAMES order for scipy.integrate.solve_ivp."""
    del t
    M1, M2, IL1, IL12, IL10, IL4, Aq, Ap, D = np.maximum(y, 0.0)
    p = params
    r = rates(y, p)
    dM1 = r["Rm1"] * r["mr"] * (1 + p.p_ros * p.ROS) - (r["Rms"] + p.mu_M1) * M1
    dM2 = r["Rm2"] * r["mr"] * max(0.0, 1 - p.i_ros * p.ROS) + r["Rms"] * M1 - p.mu_M2 * M2
    dIL1 = p.k_pn1 * r["Rp"] - p.mu_n1 * IL1
    dIL12 = p.k_pn12 * r["Rp"] - p.mu_n12 * IL12
    dIL10 = p.k_tn10 * r["Rt"] + p.k_pn10 * r["Ra"] - p.mu_n10 * IL10
    dIL4 = p.k_tn4 * r["Rt"] + p.k_pn4 * r["Ra"] - p.mu_n4 * IL4
    aq_to_ap = p.k_Aq2Ap * Aq * (p.k_Ap_M1 * M1 + p.k_Ap_D * D)
    ap_to_aq = p.k_Ap2Aq * Ap * M2
    dAq = p.s_Aq + p.k_Aq_M2 * M2 + ap_to_aq - aq_to_ap - p.mu_Aq * Aq
    dAp = aq_to_ap - ap_to_aq - p.mu_Ap * Ap
    damage = (p.alpha_n12 * IL12 + p.alpha_n1 * IL1) / (1 + (IL10 / p.a_inf2) ** 2)
    dD = damage + p.r_M1 * M1 - p.gamma_M1 * M1 * D - p.gamma_M2 * M2 * D
    dD += p.psi_Ap * Ap * D - p.psi_Aq * Aq * D
    derivatives = np.array([dM1, dM2, dIL1, dIL12, dIL10, dIL4, dAq, dAp, dD], dtype=float)
    # Retain the notebook's positivity safeguard at the boundary.
    return np.maximum(derivatives, -np.maximum(y, 0.0))


def initial_conditions():
    """Return a fresh copy of the notebook's initial state."""
    return np.array(DEFAULT_INITIAL_CONDITIONS, dtype=float)


def simulate(t_span=(0.0, 150.0), *, t_eval=None, y0=None, params=DEFAULT_PARAMETERS,
             method="RK45", rtol=1e-6, atol=1e-9):
    """Integrate the model and return scipy's ``OdeResult``."""
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 1000)
    return solve_ivp(lambda t, y: rhs(t, y, params), t_span,
                     initial_conditions() if y0 is None else np.asarray(y0, dtype=float),
                     t_eval=t_eval, method=method, rtol=rtol, atol=atol)
