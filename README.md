# Glial Crosstalk and Neuroinflammation: A Qualitative ODE Model

## Overview

This repository contains a qualitative, dimensionless model of neuroinflammatory dynamics following acute brain injury or insult. The model captures macroscopic interactions between microglial polarization states and astrocytic phenotypes in response to reactive oxygen species (ROS), with feedback mechanisms controlling neuronal damage and tissue repair.

Rather than attempting to capture molecular-scale details, this model operates at the population level: tracking the balance between pro-inflammatory M1 microglia (amplifying damage) and reparative M2 microglia (promoting recovery), alongside the recruitment and phenotypic switching of astrocytes between quiescent (Aq) and reactive (Ap) states. The interplay of these glial populations, their cytokine outputs, and ROS-dependent feedback determines whether the system resolves inflammation or transitions to chronic neurodegeneration.

## Biological Mechanisms

### Neuroinflammation and Glial Activation

Following CNS injury, resident microglia rapidly respond to danger signals, becoming activated and polarizing into distinct phenotypes. **M1 microglia** (classical/pro-inflammatory activation) produce pro-inflammatory cytokines IL-1β and IL-12, amplify ROS production, and drive neuronal damage. In contrast, **M2 microglia** (alternative/reparative activation) produce anti-inflammatory cytokines IL-4 and IL-10, which suppress M1 activity and promote tissue repair.

This polarization is not fixed: M1 microglia can be converted to M2 phenotype through exposure to anti-inflammatory signals (IL-4, IL-10), and this switch is a critical point of intervention in many therapeutic strategies.

### Astrocytic Phenotype Switching

Astrocytes exist in dynamic equilibrium between a quiescent state (Aq) and a reactive/proliferating state (Ap). **Quiescent astrocytes** maintain the resting brain environment, providing metabolic support and promoting M2 microglial recruitment. **Reactive astrocytes** respond to damage signals (M1 cytokines, ROS, physical disruption), proliferate, and can amplify inflammation by recruiting additional M1 microglia—a positive feedback loop that can drive pathology if unchecked.

This phenotypic plasticity means that even in the same injury context, the balance of astrocytic states shapes whether inflammation resolves or persists.

### ROS-Dependent Polarization

A key feature of this model is the **ROS-dependent control of microglial polarization**. Elevated ROS biases microglia toward M1 activation and away from M2 development. This creates a critical bifurcation: 
- **Low ROS scenarios** allow M2 microglia to accumulate and resolve inflammation.
- **High ROS scenarios** lock microglia into M1 phenotype, perpetuating damage and preventing recovery.

This mechanism captures a central observation in neuroinflammatory biology: the same injury signal, interpreted through the lens of oxidative stress, determines whether the outcome is recovery or chronic neurodegeneration.

### Neuronal Damage and Resolution

Neuronal damage (D) is driven primarily by M1 microglia and their cytokine products (IL-1β, IL-12), modulated by anti-inflammatory signals (IL-10) through Hill-function inhibition. Damage is resolved through microglial killing (M1 and M2 clearance) and astrocyte-mediated repair (promoted by Aq, suppressed by Ap). This balance defines the critical window for intervention: damage control is most effective when M2 recruitment and Aq-mediated repair can outpace M1-driven injury.

## Model Schematic

[**Insert model schematic diagram here** — showing microglial polarization states (M1 ↔ M2), astrocytic phenotypes (Aq ↔ Ap), cytokines (IL-1β, IL-12, IL-10, IL-4), ROS feedback, and damage dynamics (D).]

## Mathematical Model

The model state vector is **Y** = (M₁, M₂, IL-1β, IL-12, IL-10, IL-4, Aq, Ap, D), representing macroscopic populations of pro-inflammatory and reparative microglia, four cytokines, two astrocytic phenotypes, and accumulated neuronal damage.

### System of ODEs

$$\frac{dM_1}{dt} = R_{m1} m_r (1 + p_{\rm ROS} \cdot \text{ROS}) - (R_{ms} + \mu_{M1}) M_1$$

$$\frac{dM_2}{dt} = R_{m2} m_r \max(0, 1 - i_{\rm ROS} \cdot \text{ROS}) + R_{ms} M_1 - \mu_{M2} M_2$$

$$\frac{d[\text{IL-1}\beta]}{dt} = k_{pn1} R_p - \mu_{n1} [\text{IL-1}\beta]$$

$$\frac{d[\text{IL-12}]}{dt} = k_{pn12} R_p - \mu_{n12} [\text{IL-12}]$$

$$\frac{d[\text{IL-10}]}{dt} = k_{tn10} R_t + k_{pn10} R_a - \mu_{n10} [\text{IL-10}]$$

$$\frac{d[\text{IL-4}]}{dt} = k_{tn4} R_t + k_{pn4} R_a - \mu_{n4} [\text{IL-4}]$$

$$\frac{dA_q}{dt} = s_{A_q} + k_{A_q M_2} M_2 + (A_p \to A_q) - (A_q \to A_p) - \mu_{A_q} A_q$$

$$\frac{dA_p}{dt} = (A_q \to A_p) - (A_p \to A_q) - \mu_{A_p} A_p$$

$$\frac{dD}{dt} = \frac{\alpha_{n12}[\text{IL-12}] + \alpha_{n1}[\text{IL-1}\beta]}{1 + ([\text{IL-10}]/a_{\infty 2})^2} r_{M1} M_1 - (\gamma_{M1} M_1 + \gamma_{M2} M_2) D - (\psi_{A_p} A_p - \psi_{A_q} A_q) D$$

Where:
- **Astroglial phenotype conversion**: $A_q \to A_p = k_{A_q 2 A_p} A_q (k_{A_p M_1} M_1 + k_{A_p D} D)$ and $A_p \to A_q = k_{A_p 2 A_q} A_p M_2$
- **Cytokine production rates**: Rₚ and Rₜ capture M1-driven and M2-driven production, respectively; Rₐ represents astrocyte contribution.
- **Hill-function inhibition** (denominator in damage equation) captures the key mechanism: anti-inflammatory IL-10 suppresses damage proportional to its concentration.

### Parameter Choices and Qualitative Philosophy

All parameters are dimensionless and chosen to reflect qualitative biological relationships rather than exact molecular concentrations. Default values (including baseline ROS = 100) are stored in `src/glial_crosstalk/parameters.py`. The model is designed to be:

- **Robust to perturbations in parameter space** — the qualitative behavior (recovery vs. chronic inflammation) persists across reasonable parameter variations.
- **Interpretable at macroscopic scales** — state variables represent glial populations and accumulated damage, not molecular events.
- **Tunable for exploration** — users can modify parameters to explore how variations in microglial recruitment, astrocytic response, or ROS levels shift outcomes.

## Quick Start

### Installation

```bash
# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install package in editable mode
python -m pip install -e .
```

### Generate Figures

```bash
python scripts/generate_figures.py
```

This produces publication-quality figures including phase portraits, time-course simulations at varying ROS levels, and comparisons of damage trajectories under different intervention scenarios.

### Interactive Exploration

```bash
jupyter lab
```

Open `notebooks/model_walkthrough.ipynb` for a commented baseline simulation, bifurcation analysis as ROS varies, and a walkthrough of the model logic. The original analysis notebook (`notebooks/model_walkthrough_original.ipynb`) is preserved unchanged for reference.

## Repository Structure

```
glial-crosstalk-ode-model/
├── src/glial_crosstalk/
│   ├── core.py              # ODE system definition
│   ├── rates.py             # Hill functions and auxiliary dynamics
│   └── parameters.py        # Default parameter values
├── scripts/
│   └── generate_figures.py  # Publication figure generation
├── notebooks/
│   ├── model_walkthrough.ipynb        # Annotated analysis
│   └── model_walkthrough_original.ipynb # Original reference
└── README.md
```

## References and Future Work

This model builds on classical frameworks of glial activation and extends them with explicit ROS-dependent feedback. Future directions include:
- Validation against multi-photon imaging data of glial dynamics in vivo.
- Extension to spatial PDEs to capture local vs. systemic inflammatory gradients.
- Integration of additional damage pathways (excitotoxicity, metabolic failure).

---

**Questions or suggestions?** Open an issue or contact the repository maintainer.
