#!/usr/bin/env python3
"""Standalone Euler simulation extracted from the notebook example.

Run from the repository root with:
    python scripts/euler_microglial_polarization.py --ros 2.0 --steps 2000 --dt 0.01
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def hill(value: float, threshold: float, exponent: float) -> float:
    value = max(float(value), 0.0)
    threshold = max(float(threshold), 0.0)
    numerator = value ** exponent
    denominator = threshold ** exponent + numerator
    return 0.0 if denominator == 0.0 else numerator / denominator


def calculate_rm1(
    IL1,
    IL12,
    IL4,
    IL10,
    Ap,
    Aq,
    *,
    b_n,
    x_n,
    a_inf1,
    k_n1,
    k_n12,
    beta_Ap_M1,
    beta_Aq_inhib_M1,
):
    rm1 = hill(k_n1 * IL1 + k_n12 * IL12, b_n, x_n)
    rm1 *= (1.0 + beta_Ap_M1 * Ap) / (1.0 + beta_Aq_inhib_M1 * Aq)
    rm1 /= 1.0 + ((IL10 + IL4) / a_inf1) ** 2
    return rm1


def calculate_rm2(IL4, IL10, Aq, *, y_n, z_n, k_n4, k_n10, beta_Aq_M2):
    rm2 = hill(k_n4 * IL4 + k_n10 * IL10, y_n, z_n)
    rm2 *= 1.0 + beta_Aq_M2 * Aq
    return rm2


def calculate_rms(IL4, IL10, *, tau_n4, tau_n10, g_n, m_n):
    return hill(tau_n4 * IL4 + tau_n10 * IL10, m_n, g_n)


def calculate_mr(Rm1, Rm2, *, s_mr, mu_mr):
    return s_mr / (Rm1 + Rm2 + mu_mr)


def calculate_dM1_dt(Rm1, mr, Rms, M1, r_M1, gamma_M1, alpha_n1, mu_M1, ROS, p_ros):
    production = Rm1 * mr * (1.0 + p_ros * ROS)
    removal_conversion = Rms * M1 + mu_M1 * M1
    return production - removal_conversion


def calculate_dM2_dt(Rm2, mr, Rms, M1, M2, gamma_M2, mu_M2, ROS, i_ros):
    production_conversion = Rm2 * mr * max(0.0, 1.0 - i_ros * ROS) + Rms * M1
    degradation = mu_M2 * M2
    return production_conversion - degradation


def run_euler_simulation(
    *,
    ros=2.0,
    dt=0.01,
    n_steps=2000,
    il1=0.10,
    il12=0.10,
    il10=0.10,
    il4=0.10,
    aq=0.10,
    ap=0.10,
    d=0.10,
    m1=0.10,
    m2=0.10,
):
    """Integrate the simple notebook Euler model for M1/M2 dynamics."""
    params = {
        "b_n": 0.8682932268535773,
        "x_n": 2.47585316397624,
        "a_inf1": 2.944948653936038,
        "k_n1": 1.23121468628907,
        "k_n12": 1.3444424885940314,
        "beta_Ap_M1": 0.4,
        "beta_Aq_inhib_M1": 0.3,
        "y_n": 0.9996433532137823,
        "z_n": 2.3786552539888537,
        "k_n4": 1.272493506699988,
        "k_n10": 0.8577380244036067,
        "beta_Aq_M2": 0.5,
        "tau_n4": 1.1705406202322937,
        "tau_n10": 0.9048457485245999,
        "g_n": 2.464228895378905,
        "m_n": 0.563495322717254,
        "s_mr": 1.0960359353301747,
        "mu_mr": 0.36232256985896183,
        "mu_M1": 0.01,
        "mu_M2": 0.01,
        "p_ros": 1.0,
        "i_ros": 10.0,
    }

    time = np.arange(n_steps, dtype=float) * dt
    m1_trace = np.zeros(n_steps, dtype=float)
    m2_trace = np.zeros(n_steps, dtype=float)
    m1_trace[0] = m1
    m2_trace[0] = m2

    current_m1 = float(m1)
    current_m2 = float(m2)
    current_il1 = float(il1)
    current_il12 = float(il12)
    current_il10 = float(il10)
    current_il4 = float(il4)
    current_aq = float(aq)
    current_ap = float(ap)

    for step in range(1, n_steps):
        rm1 = calculate_rm1(
            current_il1,
            current_il12,
            current_il4,
            current_il10,
            current_ap,
            current_aq,
            b_n=params["b_n"],
            x_n=params["x_n"],
            a_inf1=params["a_inf1"],
            k_n1=params["k_n1"],
            k_n12=params["k_n12"],
            beta_Ap_M1=params["beta_Ap_M1"],
            beta_Aq_inhib_M1=params["beta_Aq_inhib_M1"],
        )
        rm2 = calculate_rm2(
            current_il4,
            current_il10,
            current_aq,
            y_n=params["y_n"],
            z_n=params["z_n"],
            k_n4=params["k_n4"],
            k_n10=params["k_n10"],
            beta_Aq_M2=params["beta_Aq_M2"],
        )
        rms = calculate_rms(
            current_il4,
            current_il10,
            tau_n4=params["tau_n4"],
            tau_n10=params["tau_n10"],
            g_n=params["g_n"],
            m_n=params["m_n"],
        )
        mr = calculate_mr(rm1, rm2, s_mr=params["s_mr"], mu_mr=params["mu_mr"])

        d_m1 = calculate_dM1_dt(
            rm1,
            mr,
            rms,
            current_m1,
            r_M1=0.26862120846959575,
            gamma_M1=0.38999592743353073,
            alpha_n1=1.1915264432034482,
            mu_M1=params["mu_M1"],
            ROS=ros,
            p_ros=params["p_ros"],
        )
        d_m2 = calculate_dM2_dt(
            rm2,
            mr,
            rms,
            current_m1,
            current_m2,
            gamma_M2=0.5308181970140535,
            mu_M2=params["mu_M2"],
            ROS=ros,
            i_ros=params["i_ros"],
        )

        current_m1 += dt * d_m1
        current_m2 += dt * d_m2
        current_m1 = max(current_m1, 0.0)
        current_m2 = max(current_m2, 0.0)

        m1_trace[step] = current_m1
        m2_trace[step] = current_m2

    return time, m1_trace, m2_trace


def main():
    parser = argparse.ArgumentParser(description="Euler simulation of M1/M2 microglial polarization.")
    parser.add_argument("--ros", type=float, default=2.0, help="ROS input used for the example model.")
    parser.add_argument("--dt", type=float, default=0.01, help="Euler time step in arbitrary units.")
    parser.add_argument("--steps", type=int, default=2000, help="Number of Euler steps.")
    parser.add_argument("--output", type=str, default=None, help="Optional output PNG path.")
    parser.add_argument("--no-show", action="store_true", help="Do not call plt.show().")
    args = parser.parse_args()

    time, m1_trace, m2_trace = run_euler_simulation(ros=args.ros, dt=args.dt, n_steps=args.steps)

    plt.figure(figsize=(8, 4))
    plt.plot(time, m1_trace, label="M1 (pro-inflammatory)")
    plt.plot(time, m2_trace, label="M2 (anti-inflammatory)")
    plt.xlabel("Time (arbitrary units)")
    plt.ylabel("Concentration (arbitrary units)")
    plt.title(f"Example M1/M2 dynamics with ROS={args.ros}")
    plt.legend()
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=200)
    if not args.no_show:
        plt.show()
    plt.close()

    print(f"ROS={args.ros} dt={args.dt} steps={args.steps}")
    print(f"Final M1={m1_trace[-1]:.6f}, Final M2={m2_trace[-1]:.6f}")


if __name__ == "__main__":
    main()
