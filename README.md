# The MOND interpolation function and flat rotation curves are mutually exclusive for minimally coupled k-essence

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22953268.svg)](https://doi.org/10.5281/zenodo.22953268)

**J. P. Figueroa Torres** · Independent Researcher, Guadalajara, Jalisco, Mexico · ORCID [0009-0005-5297-8777](https://orcid.org/0009-0005-5297-8777)

## Result

For a k-essence scalar with kinetic function `F(X)` and **minimally coupled** matter (the scalar acts on test particles only through the back-reaction of its stress-energy tensor), in the static spherical exterior:

- For `F = A X^n`: `v² ∝ r^(-2/(2n-1))` (for `n > 3/2`). Flat curves appear only as `n → ∞`, an isothermal dark-matter branch.
- The deep-MOND limit `μ(x) → x` requires `n = 3/2`, which gives `v² ∝ ln r / r` — not flat.
- **Corollary:** the exclusion holds for any `F` whose infrared branch is `X^(3/2)`.
- Exact flatness for arbitrary `F` fixes `F = A(√X − λ)²` uniquely, which carries a ghost.

The flat-curve configurations of Armendariz-Picon & Lim (JCAP 2005) are analysed family by family; none reaches the deep-MOND regime, and their generalised Chaplygin family reproduces the scaling law above exactly.

## Files

| File | Content |
|---|---|
| `nogo_kessence.pdf` | The paper (7 pp., 3 figures) |
| `DOS_TEOREMAS_2026-09-24.py` | Self-contained symbolic verification — **31/31 checks pass** |
| `figuras_nogo.py` | Generates the three figures of the paper |

## Reproduce

```bash
pip install sympy numpy matplotlib
python3 DOS_TEOREMAS_2026-09-24.py   # prints "31/31 PASS"
python3 figuras_nogo.py              # writes fig1_exponente.pdf, fig2_curvas.pdf, fig3_exacta.pdf
```

## Cite

J. P. Figueroa Torres, *The MOND interpolation function and flat rotation curves are mutually exclusive for minimally coupled k-essence*, Zenodo preprint (2026), [doi:10.5281/zenodo.22953268](https://doi.org/10.5281/zenodo.22953268).

## License

CC-BY-4.0.
