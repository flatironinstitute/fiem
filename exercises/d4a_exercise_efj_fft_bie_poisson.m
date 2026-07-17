%D4A_EXERCISE_EFJ_FFT_BIE_POISSON Day 4 afternoon exercise.
%
% A two-dimensional Poisson solver: EFJ extension, FFT, and BIE correction
%
% Implement the central workflow for the first two-dimensional example in
% the EFJ function-extension paper:
%
%   source f on Omega -> EFJ extension f_e on a rectangle
%   -> FFT particular solution u_p -> DLP correction u_h
%   -> u = u_p + u_h.
%
% The point of the exercise is that a smooth, windowed extension vanishing at
% the rectangle boundary has a smooth periodic continuation.  A standard FFT
% Poisson solve is therefore sufficient for the particular solution; no
% truncated-kernel method is used here.
%
% AI tools may be used to write MATLAB code.  You are responsible for checking
% the zero Fourier mode, the boundary data for the harmonic correction, and
% every reported error.
%
% Required software and mathematical sources
% ------------------------------------------
% 1. chunkie, for the boundary integral equation and its corrected evaluation:
%       https://github.com/fastalgorithms/chunkie
%
% 2. FINUFFT MATLAB interface:
%       https://finufft.readthedocs.io/en/stable/matlab.html
%
% 3. The Day 1 and Day 4 course notes, plus the paper cited below.
%
% Course starter functions included with these exercises:
%
%   d4a_helper_efj_disk_extension.m
%   d4a_helper_periodic_poisson_fft.m
%
% Reference:
%   C. L. Epstein, F. Fryklund, and S. Jiang, "An accurate and efficient
%   scheme for function extension on smooth domains," SIAM Journal on
%   Numerical Analysis 63(4) (2025), 1427--1453.

%% 1. Manufactured Poisson problem on the unit disk
%
% Let Omega = {x in R^2 : |x| < 1}.  Use the paper's first two-dimensional
% example:
%
%   Delta u = f                         in Omega,
%   u = g                               on boundary(Omega),
%
%   f(x1,x2) = -100*sin(16*pi*x1)*sin(16*pi*x2),
%
%   g(x1,x2) = 50*sin(16*pi*x1)*sin(16*pi*x2)/(16*pi)^2.
%
% The exact solution is u_ex = g throughout Omega.  Check directly that
% Delta u_ex = f before coding.  Do not use the current `rhs.example1.m`
% helper without inspecting it: define f and g above in your script so that
% the numerical data agree with the published benchmark.
%
% Use a chunkie circle with a Gauss-Legendre panel discretization.  Retain the
% same curve and boundary nodes throughout the grid-refinement study.  Use a
% modest number of panels and grid resolutions for the worksheet; the long
% 25-resolution study in the paper is not required.

%% 2. EFJ extension of f to an enclosing square
%
% On a uniform grid in B = [-L,L]^2, label the nodes in Omega and set f to
% zero at all other nodes before extension.  Read `d4a_helper_efj_disk_extension.m`,
% then use it as the starter implementation of the EFJ formula.  For an
% exterior point x = y+d*n_y, where y=x/|x| is its boundary footpoint, n_y=y,
% and d=|x|-1, it sets
%
%   f_e(x) = Phi(d) * sum_{j=0}^n w_j f(y-t_j*d*n_y),
%
% while f_e(x)=f(x) in Omega and f_e(x)=0 when d >= r1.  It computes
%
%   t_j = (a/2)*(1-cos(j*pi/n)),          w_j = ell_j(-1),
%
% as in the Day 4 morning notes, with ell_j the Lagrange basis polynomial for
% the nodes {t_j}.  Do not solve a Vandermonde system numerically.
%
% For the initial calculation, use the paper parameters
%
%   norder = 8;
%   a = 1;
%   r0 = 1e-6;
%   r1 = 32*h,
%
% and a smooth cutoff Phi that is one near d=0 and zero for d >= r1.  The
% starter uses a self-contained C-infinity cutoff with these properties,
% rather than the paper's PSWF window.  Here h is the actual uniform-grid
% spacing.  The window must make f_e negligible at the boundary of B; this is
% what permits a smooth periodic continuation for the FFT step.  Do not
% present this course calculation as a numerical reproduction of the paper's
% PSWF-window results.
%
% Make one plot of f in Omega and its extension f_e in B.  Also report
%
%   max_{x on boundary(B)} |f_e(x)|
%
% using the outer grid rows and columns.  If this value is not small compared
% with max |f_e|, stop and repair the extension/window setup before using the
% FFT.

%% 3. Periodic FFT particular solution and the zero Fourier mode
%
% Let f0 be the discrete mean of f_e on B.  Split
%
%   f_e = (f_e-f0) + f0.
%
% The mean-zero part has the periodic Fourier solution
%
%   u_per_hat(k) = -f_e_hat(k)/|k|^2,       k != 0,
%   u_per_hat(0) = 0.
%
% The missing mean is handled by the explicit particular solution
%
%   u_0(x1,x2) = f0*x1^2/2,
%
% since Delta u_0 = f0.  Thus use
%
%   u_p = u_per + u_0.
%
% Read and use `d4a_helper_periodic_poisson_fft.m`, which implements this centered
% `fft2`/`ifft2` calculation.  Apply the SAME quadratic term when evaluating
% u_p on boundary nodes.  Do not silently discard f0.  Save the Fourier
% coefficients of u_per; the quadratic term is not periodic and is not
% represented by those coefficients.

%% 4. Accurate boundary values of the particular solution
%
% The DLP correction needs u_p on the Gauss-Legendre boundary nodes, not only
% at Cartesian grid nodes.  Evaluate the periodic Fourier series u_per at
% those nodes with FINUFFT type 2, then add u_0(x1,x2) there explicitly.
%
% Let
%
%   g_h = g - u_p|_boundary(Omega).
%
% Do not replace this calculation by nearest-grid interpolation.  The boundary
% data determine the BIE correction, so interpolation error there directly
% contaminates the final solution.

%% 5. Harmonic DLP correction
%
% Solve
%
%   Delta u_h = 0                       in Omega,
%   u_h = g_h                           on boundary(Omega)
%
% with chunkie.  With the outward normal and the Day 1 interior trace
% convention, derive and solve
%
%   (-1/2*I + D) mu = g_h.
%
% Evaluate u_h at the interior grid nodes using chunkie's native corrected
% evaluation for targets close to the boundary.
%
% Form
%
%   u_num = u_p + u_h
%
% at all grid points in Omega.  Record the linear-system residual returned or
% reconstructed from the DLP solve, together with
%
%   ||u_num-u_ex||_2 / ||u_ex||_2,
%   ||u_num-u_ex||_inf / ||u_ex||_inf.
%
% Do not substitute direct smooth quadrature for targets close to the boundary.

%% 6. Resolution and extension-order studies
%
% (a) With norder=8, run at least four increasing uniform-grid resolutions.
%     Keep r1=32*h and all other geometric parameters fixed.  Tabulate h, the
%     grid size, the extension boundary check from Section 2, both solution
%     errors, and the BIE residual.
%
% (b) At one grid resolution that is not dominated by coarse-grid error, run
%     at least three extension orders, including norder=8.  Use the EFJ nodes
%     and weights supplied by the reference implementation for each order.
%     Report the same checks and plot the extended sources together if their
%     differences are visible.
%
% Plot error versus h for part (a).  Do not fit or claim an asymptotic rate
% from data that already show saturation or are clearly pre-asymptotic.

%% 7. Submit
%
% Submit one script or live script, the extension plot, the solution-error
% plot, and the two tables.  In a short explanation, answer:
%
%   * Why does the windowed EFJ extension allow the regular periodic FFT?
%   * Why is the quadratic u_0 term required when f0 is nonzero?
%   * Why is u_p|_boundary needed before solving the DLP equation?
%   * Which portions of the calculation depend on the volume grid, and which
%     depend on the boundary discretization?
%
% Record the MATLAB, FINUFFT, and chunkie revisions.
% Do not report hardware timings as general algorithmic performance results.
