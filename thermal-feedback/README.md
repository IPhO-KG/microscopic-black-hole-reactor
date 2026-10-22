# Paper II — Thermal Feedback in a Microscopic Black Hole Reactor

**Paper:** [Read the PDF](main.pdf) · Revised 22 October 2026  
**Author:** Emir Anaraly uulu, Bishkek, Kyrgyz Republic  
**Builds on:** [Paper I — fixed-medium model](../paper/)

Paper I showed that the fixed-medium Hawking–Bondi balance is unstable. Here the medium can heat up: Hawking power deposited in a reservoir raises an effective temperature $T$, which raises the sound speed and lowers the Bondi coefficient, $A(T)\propto T^{-3/2}$. In dimensionless form,

$$\frac{dx}{ds}=\theta^{-3/2}x^2-x^{-2},\qquad \frac{d\theta}{ds}=\eta x^{-2}-\zeta(\theta-1),\qquad \chi=\eta/\zeta .$$

**Results**

- There is exactly one positive equilibrium, hotter and more massive than in Paper I.
- Its Jacobian determinant is negative for every $\eta,\zeta>0$, so it is a **saddle**. The same holds for any one-zone model in which heating falls with mass, accretion falls with temperature and cooling rises with temperature.
- The physical reason is the sign of the loop: a Hawking-powered hole radiates more as it loses mass, so heating suppresses accretion exactly when it is needed.
- Timescale estimates for two hydrogen plasmas show that gas crosses the Bondi region faster than it heats or cools, so neither supports a one-zone reservoir.

## Files

| Path | Purpose |
| --- | --- |
| [`main.pdf`](main.pdf) | Manuscript |
| [`main.tex`](main.tex) | LaTeX source (the feedback-loop diagram is drawn in TikZ) |
| [`reproduce.py`](reproduce.py) | Equilibria, eigenvalues, Tables 1–2, parameter scans and Figures 2–5 |
| [`results.json`](results.json) | Table values and numerical checks |
| [`figures/`](figures/) | Vector PDF figures |

## Reproduce

From the repository root, after installing `requirements.txt`:

```bash
cd thermal-feedback
python reproduce.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The script prints Tables 1 and 2 and checks that, over a 41×41 grid with $10^{-2}\le\chi,\zeta\le10^{2}$, the determinant is negative, $\lambda_+>a>0$ and $\lambda_-<-\zeta$ at every point.

## Scope

The result is a theorem about the stated reduced model. It does not describe a realizable reactor or a hydrodynamic accretion flow.
