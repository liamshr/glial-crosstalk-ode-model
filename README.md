# Qualitative ODE Model of ROS Driven Glial crosstalk and Neuroinflammation

## Overview

This repository contains a qualitative, dimensionless model of neuroinflammatory dynamics following acute brain injury or inflammatory insult. The model captures macroscopic cellular interactions between microglial and astrocytic inflammatory polarization states in response to reactive oxygen species (ROS).

This model, rather than attempting to capture molecular-scale details of the feedback mechanisms, operates at the population level, tracking the balance between pro-inflammatory M1 microglia (amplifying damage) and anti-inflammatory M2 microglia (promoting recovery), within a relative threshold of reasonable pro-inflammatory reactions, alongside the recruitment of astrocytes between quiescent (Aq) and reactive pro-inflammatory (Ap) states. The interplay of these glial populations, their cytokine outputs, and ROS-dependent feedback determines whether the system resolves to normalized inflammation levels or transitions to chronic neurodegeneration.

## Biological Mechanisms

### Neuroinflammation and Glial Activation

Following injury in the central nervous system, microglia rapidly respond to signals of cellular threat, activating and polarizing into distinct phenotypes. **M1 microglia** produce pro-inflammatory cytokines such as IL-1β and IL-12, which amplify ROS production, and drive neuronal damage. Contrastingly, **M2 microglia** produce anti-inflammatory cytokines IL-4 and IL-10, which suppress M1 activity and promote tissue repair.

This polarization is dynamic as M1 microglia can be converted to the M2 phenotype through exposure to anti-inflammatory signals (IL-4 and IL-10 cytokines), and vice versa.

### Astrocytic Phenotype Switching

Astrocytes exist dynamically between a quiescent state (Aq) and a reactive/proliferating state (Ap), more in line with the pro-inflammatory state of microglia. **Quiescent astrocytes** maintain the resting brain environment, providing metabolic support for neurons and promoting M2 microglial polarization. **Reactive astrocytes** respond to damage signals (M1 cytokines, ROS, physical disruption), proliferate, and can amplify inflammation by recruiting additional M1 microglia. This key interaction between glial cells, although built on basic support mechanisms, can create a positive feedback loop that can drive disease progression.

Because astrocytes are so adaptable, the ratio of astrocytic functional stress determines whether tissue inflammation resolves or becomes chronic after an injury.

### ROS-Dependent Polarization

A key feature of this model is the **ROS-dependent control of microglial polarization**. Elevated ROS biases microglia toward M1 activation and away from M2 development. This creates a critical bifurcation in the system: 
- **Low ROS scenarios** allow M2 microglia to accumulate and resolve inflammation.
- **High ROS scenarios** polarize microglia into M1 phenotype, sustaining damage and preventing recovery.

This mechanism reflects a central observation in neuroinflammatory biology: the same injury signal can lead to either recovery or chronic neurodegeneration, depending on how it is shaped by oxidative stress (specific in this model to the initialized level).

### Neuronal Damage

Neuronal damage (D) is driven primarily by M1 microglia and their cytokine products (IL-1β, IL-12), while being moderated by anti-inflammatory signals (IL-10) through Hill-function inhibition. Damage is resolved through microglial-mediated clearance by both M1 and M2 microglia, along with astrocyte-mediated repair (promoted by Aq and suppressed by Ap). This balance defines a critical window for therapeutic intervention: damage control is most effective when M2 recruitment and Aq-mediated repair are able to offset M1-driven injury.

## Model Schematic

[**Model schematic diagram**]

## Mathematical Model

The model state vector is $Y = [M_1, M_2, \text{IL-1}\beta, \text{IL-12}, \text{IL-10}, \text{IL-4}, A_q, A_p, D]^T$, representing the macroscopic populations of:
* **Pro- and anti-inflammatory microglia**: $M_1$ and $M_2$
* **Pro- and anti-inflammatory cytokines**: $\text{IL-1}\beta$, $\text{IL-12}$, $\text{IL-10}$, and $\text{IL-4}$
* **Pro- and anti-inflammatory astrocytic phenotypes**: $A_p$ and $A_q$
* **Accumulated neuronal damage**: $D$

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

This produces high-resolution PNG and PDF figures for model trajectories at constant
ROS inputs of 0.00, 0.01, 31.62, and 100000.00, the ROS sensitivity analysis,
and the standalone Euler microglial polarization simulation at the same four ROS
inputs. The generation script does not produce phase portraits.

### Interactive Exploration

```bash
jupyter lab
```

Open `notebooks/model_walkthrough.ipynb` for a commented baseline simulation, bifurcation analysis as ROS varies, and a walkthrough of the model logic. The original analysis notebook (`notebooks/model_walkthrough_original.ipynb`) is saved for reference.

## Repository Structure

```
glial-crosstalk-ode-model/
├── src/
│   └── glial_crosstalk/
│       ├── __init__.py
│       ├── core.py
│       ├── parameters.py
│       └── plotting.py
├── scripts/
│   ├── euler_microglial_polarization.py
│   ├── generate_figures.py
│   └── sensitivity_analysis.py
├── notebooks/
│   ├── model_walkthrough.ipynb
│   └── model_walkthrough_original.ipynb
├── figures/                 # generated trajectory, sensitivity, and Euler figures
├── README.md
└── pyproject.toml
```

## Model Origin and Attribution

The mathematical framework implemented in this repository is based primarily on the models developed by **Vaughan et al. (2018)** and **Puri et al. (2010)**. The present implementation adapts and extends these formulations to investigate how ROS-dependent microglial polarization and feedback between microglial and astrocytic phenotypes influence the progression of neuroinflammation and neurodegeneration following acute brain injury or inflammatory insult.

While the primary equations and underlying biological mechanisms are derived from these published models, several terms and parameterizations have been modified to reflect the goals and assumptions of the present model. In particular, the model extends the underlying framework by incorporating astrocyte–microglia crosstalk and associated feedback mechanisms, as well as the direct incorporation of ROS as a discrete independent input. The resulting system is presented as an independent implementation and extension of the referenced frameworks.

## Parameterization

Parameters from the underlying microglia based equations in the neuroinflammatory model were selected using the ranges and biological constraints reported by Vaughan et al. (2018). Where parameters were provided as ranges, the midpoint was used as the baseline value to provide a reproducible nominal parameter set. Parameters introduced for astrocyte-mediated interactions were assigned temporary values based on qualitative biological relationships and the intended relative timescales of the modeled processes. The overall simulations were robust and for a large amount of parameters, only changed when parameters were shifted drastically. These parameters are treated as modeling assumptions rather than experimentally validated estimates. Parameter sensitivity for the modeling goal was subsequently evaluated to assess their influence on model behavior.

## Future Research and Alternate Directions

This model builds on classical frameworks of glial activation and extends them with explicit ROS-dependent feedback. Future directions include:
- Validation against multi-photon imaging data of glial dynamics in vivo.
- Extension to spatial PDEs to capture local vs. systemic inflammatory gradients.
- Extension to Probabilistic SDEs to account for molecular noise and cellular heterogeneity found in neuroinflammatory environments. SDEs would help model biological variability, bistability, and the random thresholds that trigger chronic glial hyperactivation.

## References

Vaughan, L. E., Ranganathan, P. R., Kumar, R. G., Wagner, A. K., & Rubin, J. E. (2018). A mathematical model of neuroinflammation in severe clinical traumatic brain injury. Journal of Neuroinflammation, 15(1). https://doi.org/10.1186/s12974-018-1384-1

Puri, I. K., & Li, L. (2010). Mathematical model for the pathogenesis of Alzheimer’s disease. PLOS ONE, 5(12), e15176. https://doi.org/10.1371/journal.pone.0015176

---

**Questions or suggestions?** Open an issue or contact the repository maintainer.
