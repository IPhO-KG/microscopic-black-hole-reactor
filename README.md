# Stability of a Microscopic Black Hole Reactor

**Paper:** [Read the current PDF](paper/main.pdf)  
**Author:** Emir Anaraly uulu, Bishkek, Kyrgyz Republic

This repository contains the manuscript, LaTeX source and code for an idealized mass-balance model. The model combines Bondi accretion with a blackbody approximation to Hawking emission:

$$\dot M=AM^2-\frac{B}{M^2},\qquad A,B>0.$$

Its only positive equilibrium is $M_{\mathrm{crit}}=(B/A)^{1/4}$. It is **unstable**: a small perturbation in mass grows rather than returning to the balance. The manuscript derives an exact implicit solution, checks numerical trajectories against it and distinguishes the model's timescales. This conclusion applies to the stated fixed-coefficient equation. The water-like parameters are a formal example; realistic emission, continuum accretion and feedback need separate treatment.

## Files

| Path | Purpose |
| --- | --- |
| [`paper/main.pdf`](paper/main.pdf) | Current eight-page manuscript |
| [`paper/main.tex`](paper/main.tex) | Editable LaTeX source |
| [`paper/reproduce.py`](paper/reproduce.py) | Recomputes constants, integrates trajectories and draws the figures |
| [`paper/results.json`](paper/results.json) | Reference values and numerical consistency checks |
| [`paper/figures/`](paper/figures/) | Three vector PDF figures |
| [`requirements.txt`](requirements.txt) | Python versions used for the supplied results |

## Reproduce

Use Python 3.11 or newer and a LaTeX installation with `pdflatex`. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cd paper
python reproduce.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

On Windows, activate the environment with `.venv\Scripts\activate` before continuing. The script writes `results.json` and regenerates `figures/figure1.pdf` through `figure3.pdf`; `pdflatex` writes `main.pdf`. It checks its numerical integrations against the analytic implicit solution.

## Scope

The blackbody lifetime and water-like Bondi parameters are illustrative outputs of the stated model, not predictions for a realizable energy system. The source PDF and code are provided so the derivations and plots can be inspected and reproduced.
