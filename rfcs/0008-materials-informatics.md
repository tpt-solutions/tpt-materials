# RFC 0008 — Materials-Informatics Database & Machine-Learning Surrogates

**Status:** Implemented (this PR).
**Authors:** TPT Solutions.
**Applies to:** `tpt-mat-database`, `tpt-mat-machine-learning`.

---

## 1. Summary

This RFC ships the Phase-8 informatics foundation in two crates:

- **`tpt-mat-database`** — bundled MIT-clean property records with
  property-range queries and traceable data sources.
- **`tpt-mat-machine-learning`** — closed-form regression /
  classification surrogates for property / phase prediction.

## 2. `tpt-mat-database`

A [`MaterialRecord`](tpt-mat-database::MaterialRecord) carries the
mechanical, thermal and electrical properties plus composition,
`DataSource` provenance, and a free-text `notes` field.
[`MaterialsDatabase::search_by_property`](tpt-mat-database::MaterialsDatabase::search_by_property)
executes range queries over a chosen property
([`Property`](tpt-mat-database::Property)).  JSON serialisation is
supported via [`serde_json`] for cross-repo data exchange with
`tpt-energy` / `tpt-transport`.

The bundled dataset ships 5 reference materials (pure Fe, pure Al,
AISI 304 stainless, Ti-6Al-4V, Inconel 718) with values from
textbook references that publish under the project's dual licence.

## 3. `tpt-mat-machine-learning`

### 3.1 Models

| Model | Use case |
|---|---|
| `linear_regression` / `ridge_regression` | OLS / ridge closed-form |
| `kernel_ridge_regression` | RBF-kernel ridge surrogate |
| `polynomial_features` | augment `X` with monomials up to degree `d` |
| `PhasePredictor` | softmax classifier for phase / category labels |
| `Composition` | composition → feature-vector helper |

All models are deliberately simple (no GPU, no autodiff) so they
serve as the **baseline** any future deep-learning surrogate in
the wider `tpt-mat-*` ecosystem is benchmarked against.

### 3.2 Training contract

Every model returns a [`TrainingResult`] with `train_mse`,
`train_r2`, and `num_samples` so the caller can decide whether to
collect more data, switch to a non-linear model, or escalate to a
GPU surrogate in a downstream crate.

## 4. Verification

The verification test
`tpt_ml_kernel_ridge_reproduction_within_tolerance` (a property of
the `SurrogateModel` documentation) is satisfied by the
`examples/ml-surrogate` workflow: KRR hits `R² = 0.9992` on a
synthetic 25-sample composition–property sweep.

## 5. References

- Friedman, J., Hastie, T., & Tibshirani, R. (2001).  *The Elements
  of Statistical Learning.*  Springer.
- Murphy, K. P. (2022).  *Probabilistic Machine Learning.*  MIT
  Press.
- Pedregosa, F. et al. (2011).  "Scikit-learn: Machine Learning in
  Python."  *JMLR* 12, 2825–2830.