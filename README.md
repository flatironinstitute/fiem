# FIEM — Fast Integral Equation Methods

Materials on Fast Integral Equation Methods (FIEM): lecture notes and selected
exercises.

## Contents

- **`lecture_notes/`** — Lecture notes (LaTeX source and compiled PDFs).
- **`exercises/`** — Selected problem sets and worked solutions, including any
  accompanying code (Fortran / C / C++).

## Topics

Integral equation formulations and the fast algorithms used to solve them at
scale, including:

- Boundary and volume integral equation formulations
- Nyström and Galerkin discretizations, quadrature for singular kernels
- Fast Multipole Method (FMM) and related hierarchical / fast direct solvers
- Applications to Laplace, Helmholtz, and Stokes problems

## Contributors

- Shidong Jiang
- Manas Rachh

## Exercise Software

The exercises are MATLAB worksheets and use the following external packages.
Install the packages and make their MATLAB interfaces available on the MATLAB
path before the relevant exercise session.

- [chunkie](https://github.com/fastalgorithms/chunkie) --- MATLAB toolbox for
  integral equation discretization and solution.
- [FMM2D](https://github.com/flatironinstitute/fmm2d) --- two-dimensional fast
  multipole libraries for Laplace and Helmholtz interactions.
- [FINUFFT](https://finufft.readthedocs.io/) --- nonuniform fast Fourier
  transforms, with MATLAB interfaces.
- [FLAM](https://github.com/klho/FLAM) --- Fast Linear Algebra in MATLAB; used
  by exercises on interpolative decomposition and related compression methods.

Individual worksheets state the package requirements specific to that exercise.

## Building the notes

LaTeX sources compile with a standard TeX distribution, e.g.:

```sh
cd lecture_notes
latexmk -pdf day1_morning.tex
```

Build artifacts (LaTeX auxiliary files, compiled object files) are excluded via
`.gitignore`.
