# Stability of a Microscopic Black Hole Reactor

This repository contains the paper, source files, and numerical code for:

**Stability of a Microscopic Black Hole Reactor: Hawking Evaporation and a Fixed-Medium Bondi Model**

**Author:** Emir Anaraly uulu  
**Location:** Bishkek, Kyrgyz Republic

## Overview

This project studies whether accretion from a surrounding medium can balance Hawking evaporation for a microscopic black hole.

The model combines idealized Bondi accretion with a blackbody approximation to Hawking radiation:

\[
\dot{M} = AM^2 - \frac{B}{M^2}.
\]

The equation has one positive equilibrium mass,

\[
M_{\mathrm{crit}} = \left(\frac{B}{A}\right)^{1/4},
\]

but this equilibrium is unstable: a small increase in mass leads to further growth, while a small decrease leads to further evaporation.

The project also derives an exact implicit solution, compares it with numerical integration, and discusses the physical limitations of the model.

## Main result

Within the fixed-coefficient Hawking–Bondi model, passive mass balance is unstable.

This result applies to the stated idealized equation and is not a general impossibility result for black-hole energy systems.

## Repository contents

- `paper.pdf` — final version of the research paper
- `main.tex` — LaTeX source
- `reproduce.py` — numerical calculations and figure generation
- `results.json` — numerical results used in the paper
- `figures/` — figures used in the manuscript

## Reproducing the results

Install the required Python packages:

```bash
pip install numpy scipy matplotlib
