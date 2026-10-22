# Combined paper

**Can a Microscopic Black Hole Be Held at Constant Mass? Stability of the Hawking–Bondi balance in fixed and thermally coupled media**

[Read the PDF](main.pdf) · 22 October 2026

This manuscript presents [Paper I](../paper/) (fixed medium) and [Paper II](../thermal-feedback/) (thermal feedback) as one study, with a shared introduction, a single set of definitions and joint conclusions. It adds a short estimate of radiation pressure at the balance masses.

The figures are read directly from `../paper/figures/` and `../thermal-feedback/figures/`, so run both `reproduce.py` scripts before compiling:

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

`fixed_medium.png` and `phase_portrait.png` are raster copies of two figures for the repository front page.
