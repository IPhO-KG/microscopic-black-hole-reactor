# Paper II — Thermal Feedback in a Microscopic Black Hole Reactor

**Paper:** [Read the PDF](main.pdf) · Revised 30 September 2026  
**Author:** Emir Anaraly uulu, Bishkek, Kyrgyz Republic  
**Builds on:** [Paper I — fixed-medium model](../paper/)

Paper I showed that the fixed-medium Hawking–Bondi balance is unstable. Here the medium is allowed to heat up: Hawking power deposited into a reservoir raises an effective temperature $T$, which raises the sound speed and lowers the Bondi coefficient, $A(T)\propto T^{-3/2}$. In dimensionless form,

$$\frac{dx}{ds}=\theta^{-3/2}x^2-x^{-2},\qquad \frac{d\theta}{ds}=\eta x^{-2}-\zeta(\theta-1),\qquad \chi=\eta/\zeta .$$

**Result:** there is exactly one positive equilibrium, hotter and heavier than in Paper I, but its Jacobian determinant is negative for every $\eta,\zeta>0$ — the equilibrium is a **saddle**. The same holds for any one-zone model where heating falls with mass, accretion falls with temperature and cooling rises with temperature.

## Files

| Path | Purpose |
| --- | --- |
| [`main.pdf`](main.pdf) | Manuscript |
| [`reproduce.py`](reproduce.py) | Solves for the equilibrium, eigenvalues, Table 1, the 41×41 parameter scan and Figures 1–3 |
| [`results.json`](results.json) | Table 1 values and numerical checks |
| [`figures/`](figures/) | Vector PDF figures |

## Reproduce

From the repository root, after installing `requirements.txt`:

```bash
cd thermal-feedback
python reproduce.py
```

The script prints Table 1 and checks that, over the 41×41 grid $10^{-2}\le\chi,\zeta\le10^{2}$, the determinant is negative and $\lambda_+>a>0$ everywhere.

## Scope

This is a theorem about the stated reduced model, not a prediction for a real material or a hydrodynamic accretion flow.
