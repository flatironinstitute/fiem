%D4M_EFJ_AAA_DISK_EXTENSION Day 4 morning exercise.
%
% EFJ and AAA extensions on the unit disk
%
% Reproduce the two-dimensional disk comparison in the EFJ function-extension
% paper.  The purpose is to examine an extension as a function of two
% variables.  In particular, compare the azimuthal Fourier spectra of the
% EFJ extension and the AAA extension outside the disk.
%
% This is an extension experiment only.  Do not add the FFT particular
% solution, boundary integral correction, or a Poisson-solver timing study.
%
% AI tools may be used to write MATLAB code.  You are responsible for checking
% that every EFJ sample point lies inside the disk and that the two methods use
% identical interior data and exterior targets.
%
% Required software and mathematical sources
% ------------------------------------------
% 1. Chebfun, including its `aaa` function:
%       https://www.chebfun.org/
%
% 2. The Day 4 morning notes and the EFJ paper, for the extension formula,
%    Chebyshev nodes, and Lagrange-basis weights.
%
% References:
%   C. L. Epstein, F. Fryklund, and S. Jiang, "An accurate and efficient
%   scheme for function extension on smooth domains," SIAM Journal on
%   Numerical Analysis 63(4) (2025), 1427--1453.
%
%   Y. Nakatsukasa, O. Sete, and L. N. Trefethen, "The AAA algorithm for
%   rational approximation," SIAM Journal on Scientific Computing 40(3)
%   (2018), A1494--A1522. https://doi.org/10.1137/16M1106122

%% 1. Problem and common samples
%
% Let D be the unit disk, parameterized by the boundary point
%
%   y(theta) = (cos(theta),sin(theta)),       0 <= theta < 2*pi.
%
% The outward normal at y(theta) is y(theta) itself.  Use the first
% two-dimensional example in the EFJ paper:
%
%   f(x1,x2) = 100*sin(16*pi*x1)*sin(16*pi*x2),      (x1,x2) in D.
%
% Treat this formula as data available only in D when constructing either
% extension.  Although it has a known analytic continuation, an extension of
% a general C^n function is not required to reproduce such a continuation.
%
% Use the comparison's radial mesh
%
%   nr = 160;
%   rho = 1.2*(0:nr)'/nr;
%   rho_in  = rho(rho <= 1);
%   rho_out = rho(rho > 1);
%
% and begin with 500 equally spaced values of theta for the timing comparison.
% For the spectral plots, use the angular resolution in `testaaa2dPS.m`,
% namely ntheta = 2^12.  Use the same theta and rho arrays for EFJ and AAA.

%% 2. EFJ extension
%
% Use the eighth-order EFJ formula with
%
%   n = 8;
%   a = 2;
%
% Compute the shifted/scaled Chebyshev nodes and the weights as in the Day 4
% morning notes:
%
%   t_j = (a/2)*(1-cos(j*pi/n)),          j = 0,...,n,
%   w_j = ell_j(-1),
%
% where ell_j is the Lagrange basis polynomial for the nodes {t_j}.  Do not
% solve a Vandermonde system numerically.
%
% For an exterior target x = y(theta) + d*y(theta), with d = rho-1 > 0, use
% the formula from the paper:
%
%   E_EFJ[f](y(theta)+d*y(theta))
%       = sum_{j=0}^n w_j f((1-t_j*d)*y(theta)).
%
% This is the disk specialization of sampling at y - t_j*d*n_y.  Verify
% before evaluation that
%
%   0 <= 1-t_j*d <= 1
%
% for all samples.  With a=2 and rho <= 1.2, the smallest sample radius is
% 0.6.  Sampling at 1+t_j*d would be outside D and is not the EFJ extension
% formula.
%
% Build a matrix U_EFJ(theta,rho).  Copy the interior data into the columns
% rho <= 1 and use the formula above only for rho > 1.  Do not apply a window
% in this exercise: the comparison isolates the unwindowed extensions.

%% 3. AAA extension
%
% At each fixed theta, sample the same interior radial data
%
%   F_theta(rho) = f(rho*cos(theta),rho*sin(theta)),  0 <= rho <= 1,
%
% construct the Chebfun AAA rational approximation, and evaluate it at the
% common exterior radii:
%
%   rat = aaa(F_theta,rho_in);
%   E_AAA[f](theta,rho_out) = rat(rho_out).
%
% Assemble U_AAA(theta,rho) with the identical interior columns.  This is
% called the AAA extension.  On the disk, the normal coordinate is radial,
% but do not rename either method based on that special geometry.
%
% Time only the construction and evaluation of each extension over the 500
% angular samples.  Run each timing at least three times after a warm-up and
% report the median.  Record MATLAB, Chebfun, and reference-code revisions.
% Hardware- and version-dependent times are observations, not universal
% performance claims.

%% 4. Azimuthal smoothness diagnostic
%
% For each of the three radii
%
%   rho_star = [1.005, 1.035, 1.065],
%
% extract E(theta) from U_EFJ and U_AAA.  Compute its discrete Fourier
% coefficients in theta using the same normalization for both methods, e.g.,
%
%   Ehat = fftshift(fft(E))/ntheta;
%   k = -ntheta/2:ntheta/2-1;
%
% Make one semilogarithmic plot of abs(Ehat) versus k for each rho_star, with
% the EFJ and AAA curves together.  These are the spectra shown in the paper's
% disk comparison.  Use axes and legends large enough to reveal the tail.
%
% Also plot the extension values versus theta at rho = 1.065.  Do not infer
% smoothness from a visually smooth low-resolution curve; use the Fourier tail
% as the diagnostic.

%% 5. Report and interpretation
%
% Submit:
%
%   1. A script or live script that reproduces the common disk data, the EFJ
%      construction, and the AAA construction.
%   2. A table with n, a, nr, ntheta, the number of EFJ function evaluations
%      per exterior target, the median timing for each method, and software
%      versions.
%   3. The three azimuthal-spectrum plots and the rho=1.065 value plot.
%   4. A short interpretation addressing the questions below.
%
% Questions:
%
%   * Why does the EFJ formula guarantee C^n matching across the boundary for
%     a smooth domain, while each AAA fit is a separately selected rational
%     function?
%   * What do the azimuthal spectra test that a one-dimensional radial error
%     plot does not?
%   * In this experiment, how do the EFJ and AAA timings depend on the number
%     of evaluations used to form an exterior target value?
%
% Report the behavior shown by your data.  Do not claim that either method is
% universally faster or more accurate outside this prescribed comparison.
