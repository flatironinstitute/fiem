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

## Building the notes

LaTeX sources compile with a standard TeX distribution, e.g.:

```sh
cd lecture_notes
latexmk -pdf day1_morning.tex
```

Build artifacts (LaTeX auxiliary files, compiled object files) are excluded via
`.gitignore`.
