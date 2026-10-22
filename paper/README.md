# Paper I — Stability of a Microscopic Black Hole Reactor

**Hawking evaporation and a fixed-medium Bondi model**

**Paper:** [Read the PDF](main.pdf) · Revised 22 October 2026  
**Continued in:** [Paper II — thermal feedback](../thermal-feedback/) · **Combined version:** [full paper](../full-paper/)

The model combines Bondi accretion with a blackbody approximation to Hawking emission:

$$\dot M=AM^2-\frac{B}{M^2},\qquad A,B>0.$$

Its only positive equilibrium, $M_{\mathrm{crit}}=(B/A)^{1/4}$, is **unstable**: a small change in mass grows instead of returning to the balance. The paper derives an exact implicit solution, checks numerical trajectories against it and compares the model's timescales. The water-like parameters set numerical scales only.

## Files

| Path | Purpose |
| --- | --- |
| [`main.pdf`](main.pdf) | Eight-page manuscript |
| [`main.tex`](main.tex) | LaTeX source |
| [`reproduce.py`](reproduce.py) | Recomputes constants, integrates trajectories and draws the figures |
| [`results.json`](results.json) | Reference values and numerical checks |
| [`figures/`](figures/) | Three vector PDF figures |

## Reproduce

From the repository root, after installing `requirements.txt`:

```bash
cd paper
python reproduce.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The script writes `results.json`, regenerates the figures and checks the numerical integrations against the exact implicit solution.
