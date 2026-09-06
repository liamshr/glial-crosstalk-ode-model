# Glial crosstalk ODE model

This repository organizes the original `notebooks/model_walkthrough_original.ipynb`
into a reusable scientific Python model of neuroinflammation. It is a qualitative,
dimensionless model: activated M1 microglia amplify IL-1β/IL-12 and tissue damage,
whereas M2 microglia, IL-4/IL-10, and quiescent astroglia provide resolution and
repair. Proliferating astroglia feed back positively on M1 recruitment and can
increase damage; quiescent astroglia promote M2 recruitment and repair.

## Model

The state is \(Y=(M_1,M_2,IL1\beta,IL12,IL10,IL4,A_q,A_p,D)\). The notebook's
Hill functions are retained in `glial_crosstalk.core.rates`. The central dynamics
are:

\[
\begin{aligned}
\dot M_1 &= R_{m1}m_r(1+p_{\rm ros}ROS)-(R_{ms}+\mu_{M1})M_1,\\
\dot M_2 &= R_{m2}m_r\max(0,1-i_{\rm ros}ROS)+R_{ms}M_1-\mu_{M2}M_2,\\
\dot{IL1\beta} &= k_{pn1}R_p-\mu_{n1}IL1\beta,\quad
\dot{IL12}=k_{pn12}R_p-\mu_{n12}IL12,\\
\dot{IL10} &= k_{tn10}R_t+k_{pn10}R_a-\mu_{n10}IL10,\\
\dot{IL4} &= k_{tn4}R_t+k_{pn4}R_a-\mu_{n4}IL4,\\
\dot A_q &= s_{Aq}+k_{Aq\_M2}M_2+A_p\!\to A_q-A_q\!\to A_p-\mu_{Aq}A_q,\\
\dot A_p &= A_q\!\to A_p-A_p\!\to A_q-\mu_{Ap}A_p,\\
\dot D &= \frac{\alpha_{n12}IL12+\alpha_{n1}IL1\beta}{1+(IL10/a_{\infty2})^2}
r_{M1}M_1-(\gamma_{M1}M_1+\gamma_{M2}M_2)D
(\psi_{Ap}A_p-\psi_{Aq}A_q)D.
\end{aligned}
\]

The astroglial conversions are \(A_q\!\to A_p=k_{Aq2Ap}A_q(k_{Ap\_M1}M_1+
k_{Ap\_D}D)\) and \(A_p\!\to A_q=k_{Ap2Aq}A_pM_2\). All numerical defaults,
including the original \(ROS=100\), live in `src/glial_crosstalk/parameters.py`.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python scripts/generate_figures.py
jupyter lab
```

Open `notebooks/model_walkthrough.ipynb` for a commented baseline simulation,
phase portrait, and a comparison of microglial activation/damage trajectories
as ROS changes. The original notebook is preserved unchanged.
