# Changelog

All notable changes to this project are documented here.

## [3.0.0] – 2026-07-15

Full production release of the Model Scaling module.

### Added

**Model Scaling (`gurobi_modelanalyzer.scaling`)**

- `scale_model(model, method, ...)` — scale a Gurobi model using one of three
  iterative methods (`equilibration`, `geometric_mean`, `arithmetic_mean`) and
  return a `ScaledModel` ready for optimization.
- `read_scaling_file(path, model)` — load per-variable and per-constraint
  scaling factors from a `.scl` file and apply them to a model before calling
  `scale_model`.
- `ScaledModel.write_scaling(path, lock_factors=True)` — export the computed
  scaling factors to a `.scl` file for later reuse or warm-starting.
- `ScaledModel.getVarsUnscaled()`, `getConstrsUnscaled()`,
  `getQConstrsUnscaled()` — access solution values and violations in the
  original (unscaled) variable space.
- `ScaledModel.computeUnscVio()`, `computeUnscObj()` — compute constraint and
  bound violations and objective value in the original space.
- `ScaledModel.ColScaling`, `ScaledModel.RowScaling` — diagonal sparse scaling
  matrices.
- `scaling_factor` attribute on `ScaledVar`, `ScaledConstr`, `ScaledQConstr`
  wrapper objects.
- `power_of_2=True` option for `scale_model` — round all final scaling factors
  to the nearest power of 2 for exact floating-point representation.
- Automatic inheritance of variable attributes (`Start`, `VarHintVal`, `PStart`,
  `VarHintPri`, `BranchPriority`, `Partition`) and constraint attributes
  (`Lazy`) from the original model.
- Simplex basis (`VBasis`, `CBasis`) transferred via `.bas` file, supporting
  both user-provided warm-start bases and post-solve basis inheritance.
- `ValueError` raised for models with general constraints ((MI)NLP), which are
  not supported by the scaling module.

**`gurobi_cls` command-line tool**

- Scale and solve a Gurobi model in one step, producing a unified log (scaling
  log + Gurobi solve log + solution quality statistics) and an automatic `.scl`
  output file.
- `--input-file PATH` — load additional input files (MIP starts `.mst`,
  attribute files `.attr`, basis files `.bas`, etc.) via `model.read()`.
- `--scaling-file PATH` — supply a `.scl` file as a warm-start for the scaling
  algorithm.
- `--method`, `--scale-passes`, `--scale-conv-tol`, `--scaling-lb`,
  `--scaling-ub`, `--value-threshold`, `--scaling-time-limit`, `--power-of-2`,
  `--no-console-log` options.

### Changed

- `read_scaling_file` is now exported from the top-level
  `gurobi_modelanalyzer` namespace alongside `scale_model`.

---

## [3.0.0a1] – 2024

Alpha release of the Model Scaling module.

---

## [2.1.0]

- Bug fixes and improvements to the ill-conditioning explainer and solution
  checker modules.

## [2.0.0]

- Refactored results analyzer module.

## [1.0.1]

- Minor fixes.

## [1.0.0rc0]

- Release candidate for 1.0.

## [0.2.0b1] / [0.1.0b1] / [0.1.0b0]

- Early beta releases.
