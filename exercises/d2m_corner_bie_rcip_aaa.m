%D2M_CORNER_BIE_RCIP_AAA Day 2 morning exercise.
%
% One-corner Laplace problems: dyadic refinement, RCIP, and AAA-LS
%
% The purpose of this exercise is to compare three ways of resolving a
% corner.  Do not try to implement every detail of a production quadrature
% or RCIP code from scratch.  Start from the cited implementations, make
% the requested changes, and validate the result.
%
% AI tools may be used to write MATLAB code.  You must still identify which
% equation each program solves, run convergence tests, and explain the
% evidence in your own words.
%
% Required software and source code
% ---------------------------------
% 1. RCIP example code:
%       https://github.com/haiszhu/RCIP
%    The starting files are
%       test/test_lapCornerBVP_RCIP.m
%       test/rcip_geo.m
%       test/rcip_solve.m
%       helsing/demo4.m
%
% 2. Chebfun's AAA implementation:
%       https://www.chebfun.org/
%    The function is `aaa`.
%
% Read before coding:
%   J. Helsing, "Solving integral equations on piecewise smooth boundaries
%   using the RCIP method: a tutorial," arXiv:1207.6737v10.
%
%   S. Costa and L. N. Trefethen, "AAA-least squares rational approximation
%   and solution of Laplace problems," 9th European Congress of Mathematics,
%   2023, DOI: 10.4171/8ECM/16.
%
% Problem and common geometry
% ---------------------------
% Use the one-corner contour already implemented in `rcip_geo.m`:
%
%   z(s) = sin(pi*s) exp(i*theta*(s - 1/2)),  0 <= s <= 1.
%
% The parameter theta is the opening angle at z = 0.  Begin with the value
% used in the supplied test.  Use the same curve, boundary data, and target
% points for all three experiments.  The reference field and its source point
% in `test_lapCornerBVP_RCIP.m` define a manufactured Laplace problem; retain
% only targets where that field is harmonic.
%
% The RCIP part concerns a second-kind boundary integral equation
%
%   (I + lambda*K) rho = h,
%
% with a corner-local splitting K = K_star + K_circ.  In the compressed
% discretization the global unknown is the transformed density rho_tilde and
%
%   R = P_W^T (I_fin + lambda*K_star,fin)^(-1) P,
%
%   (I_coa + lambda*K_circ,coa*R) rho_tilde,coa = h_coa.
%
% Here P is prolongation and P_W^T is its quadrature-weighted restriction.
% Use the notation and operator convention of the RCIP tutorial and source
% code; do not substitute a different jump sign or layer representation.

%% 1. Establish the reference calculation
%
% (a) Clone the RCIP repository and run `test_lapCornerBVP_RCIP.m` unchanged.
%     Record the MATLAB version, the repository revision, and every path you
%     add.  Confirm that the error plot is produced.
%
% (b) Read the following functions and state their roles in two or three
%     sentences each:
%
%       rcip_geo:   coarse grid, dyadically refined grid, and local meshes;
%       rcip_solve: compressed coarse equation and the matrix R;
%       rcip_rhofin: backward reconstruction of the original fine density;
%       rcip_eval:  field evaluation at targets.
%
% (c) Plot the contour, the manufactured source point, and the target points
%     retained for error measurement.  State which side of the contour is
%     the physical domain for this test.  Exclude the source singularity from
%     every error norm.

%% 2. Direct dyadic refinement
%
% Use `helsing/demo4.m` as the direct-versus-RCIP template.  Its direct path
% solves a fully dyadically refined system, whereas its RCIP path retains a
% coarse global system.  Adapt the direct path to the boundary integral
% problem in Section 1, using the same panel order as the supplied code.
%
% Increase the number of dyadic levels nsub.  For each level, report:
%
%   * the number of boundary unknowns in the direct system;
%   * the maximum error at the selected targets;
%   * the boundary residual; and
%   * the wall-clock time for assembly plus solution.
%
% Stop when the target error stabilizes, or when the direct system is no
% longer practical on your machine.  This direct calculation is the reference
% against which the compressed calculation is checked; it is not expected to
% be the preferred solver at deep refinement.

%% 3. RCIP compression
%
% Run the RCIP calculation for the same sequence of nsub values.  Keep the
% global coarse mesh fixed and record:
%
%   * the coarse-system dimension;
%   * the dimension of the local blocks used in the recursion;
%   * the GMRES iteration count and residual;
%   * the target error; and
%   * the smooth functional used by the supplied example, if applicable.
%
% Verify numerically that the direct and RCIP fields agree to the requested
% accuracy at the selected targets.  Do not compare raw density values at the
% corner.  The original density can be singular there.  For a smooth
% functional, use the weight-corrected density
%
%   rho_hat,coa = R*rho_tilde,coa
%
% together with the coarse quadrature weights.  Reconstruct rho_fin only to
% plot the corner density or to evaluate a quantity that requires it.
%
% Make one plot with target error versus nsub for direct refinement and RCIP.
% Make a second plot showing the number of global unknowns versus nsub.

%% 4. AAA-LS rational solver
%
% Apply the AAA-least-squares (AAA-LS) construction of Costa and Trefethen to
% the SAME boundary samples and Dirichlet data.  This is the rational-solver
% comparison in this exercise.
%
% The required workflow is:
%
%   1. Sample the boundary and data, with adequate resolution near the corner.
%   2. Use `aaa` to construct a rational approximation of the sampled boundary
%      data, globally or locally near the corner as in the paper.
%   3. Discard AAA poles in the physical domain and retain poles in its
%      complement.
%   4. Use the retained poles in the real boundary least-squares solve for a
%      harmonic field u = Re(r).
%   5. Evaluate u at the same targets as in Sections 2--3.
%
% This is AAA-LS, not merely an AAA fit of the boundary samples.  The final
% real least-squares solve is required to obtain the harmonic field.  Plot the
% retained poles with the contour, then report the maximum boundary residual
% and maximum target error as the AAA tolerance or maximum degree is varied.
%
% Do not use the standard `laplace.m` Lightning Laplace Solver for this
% contour: that implementation accepts polygons and circular polygons, while
% the RCIP test contour is a general one-corner curve.  If you instead choose
% a polygonal benchmark, state the new geometry explicitly and use it for all
% three methods.

%% 5. Submit
%
% Submit a concise notebook or MATLAB script, the two requested convergence
% plots, one pole plot, and a table with these columns:
%
%   method, refinement level or AAA degree, global unknowns or basis size,
%   boundary residual, target error, and runtime.
%
% In a short conclusion, answer only these questions:
%
%   * Why does dyadic refinement resolve the corner but enlarge the global
%     linear system?
%   * What information does R retain after the fine corner scales are removed?
%   * What role do the exterior AAA poles play in the rational approximation?
%
% Do not claim a convergence rate from a short numerical table.  Report only
% the behavior visible in the data and identify any calculation that did not
% converge or could not be reproduced.
