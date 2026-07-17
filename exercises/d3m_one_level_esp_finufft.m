%D3M_ONE_LEVEL_ESP_FINUFFT Day 3 morning exercise.
%
% One-level Ewald summation with FINUFFT
%
% This exercise implements a small, periodic Coulomb-summation prototype.
% Its smooth contribution is evaluated as
%
%   type-1 NUFFT -> diagonal Fourier multiplier -> type-2 NUFFT,
%
% while its local contribution is evaluated by direct kernel summation.  The
% objective is to make this decomposition and its normalization explicit.
% This is an accuracy exercise, not an implementation or performance study of
% the complete ESP algorithm.
%
% AI tools may be used to write MATLAB code.  You are responsible for checking
% the Fourier signs, factors of L, the omitted zero mode, and the self term.
%
% Required software and sources
% -----------------------------
% 1. FINUFFT MATLAB interface:
%       https://finufft.readthedocs.io/en/stable/matlab.html
%
% 2. ESP companion code:
%       https://github.com/lu1and10/esp_spread
%    Clone it, add it to the MATLAB path, and run `init`.  It requires
%    Chebfun.  This exercise uses only
%
%       esp_split.m       PSWF kernel split;
%       ewald_direct.m    direct periodic reference.
%
% Read before coding:
%   J. Liang, L. Lu, A. H. Barnett, L. Greengard, and S. Jiang,
%   "Accelerating Molecular Dynamics Simulations using Fast Ewald Summation
%   with Prolates," Nature Communications 17, 73232 (2026).
%   https://doi.org/10.1038/s41467-026-73232-8
%
% Use the definitions returned by `esp_split` directly.  In particular, that
% companion code uses the 1/r convention, not the 1/(4*pi*r) convention.
% Do not mix the two normalizations.

%% 1. Periodic problem and PSWF split
%
% Work in the cubic periodic box Omega = [0,L)^3 with charges q_j at positions
% r_j.  Impose charge neutrality,
%
%   sum_j q_j = 0,
%
% and omit the k = 0 Fourier mode.  These choices fix the additive constant
% convention for the periodic potential.  The potential at particle i is
%
%   u_i = sum_{j,n}' q_j / |r_i-r_j-L*n|,
%
% where the prime omits j=i, n=0.
%
% For a reproducible test, begin with
%
%   rng(1)
%   L = 1;
%   N = 200;
%   r = L*rand(N,3);
%   q = randn(N,1);  q = q - mean(q);
%   rc = L/4;
%   split_tol = 1e-10;
%   S = esp_split(split_tol,rc);
%
% The choice rc < L/2 lets the minimum-image convention identify every pair
% that contributes to the local part.  Verify `abs(sum(q))` is near working
% precision.  Record `S.c`, `S.C0`, and `S.Kmax` in the report.
%
% The companion implementation provides the split
%
%   1/r = N^c(r) + F^c(r),
%
% with N^c supported on [0,rc].  Its Fourier multiplier is F_hat^c(|k|).
% The target formula is therefore
%
%   u_i^near = sum_{j,n}' N^c(|r_i-r_j-L*n|) q_j,
%
%   u_i^far = sum_{k in Z^3, k != 0} exp(-i*k.r_i)
%                 F_hat^c(|k|)/L^3
%                 sum_j exp(i*k.r_j) q_j,
%
%   u_i = u_i^near + u_i^far - q_i/(S.C0*rc).
%
% The last term removes the smooth self interaction included by the Fourier
% sum.  Check it against the `pot_self` calculation in `ewald_direct.m`.

%% 2. Local contribution: direct minimum-image calculation
%
% Write a short function `local_potential(r,q,L,S)` that loops over targets
% i and sources j.  For each pair, form the minimum-image displacement
%
%   d = r(i,:) - r(j,:);
%   d = d - L*round(d/L);
%
% and add `S.N(norm(d))*q(j)` only when i ~= j and norm(d) <= rc.  The first
% condition removes the physical self interaction; the second is harmless
% but makes the compact support visible in the code.
%
% Do not use this minimum-image shortcut if you change the experiment to
% rc >= L/2.  In that case, sum the periodic images required by the support,
% as `ewald_direct.m` does.
%
% Plot N^c(r) and F^c(r) for 0 < r <= L/2.  Mark rc on the plot and state
% which part has compact real-space support.

%% 3. Smooth contribution: two FINUFFTs and one diagonal multiplication
%
% Map particle coordinates to FINUFFT's 2*pi-periodic coordinates:
%
%   x = 2*pi*r(:,1)/L;
%   y = 2*pi*r(:,2)/L;
%   z = 2*pi*r(:,3)/L;
%
% Choose an even number of modes in each coordinate.  One sufficient initial
% choice is
%
%   nf = 2*ceil(S.Kmax*L/(2*pi)) + 2;
%
% because FINUFFT's default centered mode order is
%
%   -nf/2, ..., nf/2-1.
%
% Construct the mode array in that same order:
%
%   m = -nf/2:nf/2-1;
%   [m1,m2,m3] = ndgrid(m,m,m);
%   kmag = (2*pi/L)*sqrt(m1.^2 + m2.^2 + m3.^2);
%
% Use a type-1 transform with positive sign to form the structure factor
%
%   rho_hat(k) = sum_j q_j exp(i*k.r_j):
%
%   nufft_tol = 1e-12;
%   rho_hat = finufft3d1(x,y,z,q,+1,nufft_tol,nf,nf,nf);
%
% Form the diagonal multiplier and omit every mode outside the split's
% Fourier cutoff, including the zero mode:
%
%   multiplier = S.F_hat(kmag)/L^3;
%   multiplier(kmag == 0 | kmag > S.Kmax) = 0;
%   uhat = multiplier .* rho_hat;
%
% Finally use a type-2 transform with negative sign:
%
%   u_far = real(finufft3d2(x,y,z,-1,nufft_tol,uhat));
%
% Confirm from the FINUFFT documentation that your installation uses its
% default centered mode ordering.  If you change `modeord`, change the mode
% array consistently.  A mode-order mismatch gives a wrong answer even when
% both FINUFFT calls succeed.

%% 4. Reference calculation and checks
%
% Compute the reference using exactly the same split, not a different Ewald
% convention:
%
%   refopts = struct('compute_force',false,'image_shells',1, ...
%                    'kmax_factor',1.0);
%   [u_ref,comp_ref] = ewald_direct(r,q,S,L,refopts);
%
% Assemble your result as
%
%   u_num = u_near + u_far - q/(S.C0*rc);
%
% and report
%
%   rel_l2  = norm(u_num-u_ref.pot)/norm(u_ref.pot);
%   abs_linf = norm(u_num-u_ref.pot,inf);
%
% Also compare `u_near`, `u_far`, and the self term separately against
% `comp_ref.pot_near`, `comp_ref.pot_far`, and `comp_ref.pot_self`.  These
% three comparisons are required: a small cancellation in the final potential
% does not establish that the decomposition was implemented correctly.
%
% Make one scatter plot of u_num versus u_ref.pot, with the diagonal line.
% Make a second plot of the pointwise error against particle index or position.

%% 5. Resolution and NUFFT-tolerance study
%
% Keep L, r, q, rc, and S fixed.  Perform both studies below.
%
% (a) Fourier resolution.  Let nf_full be the value in Section 3.  Run at
%     least three even values ending at nf_full, including one that does not
%     represent all modes with |k| <= S.Kmax.  For each value, construct the
%     matching centered mode grid and report rel_l2 and abs_linf.  Explain the
%     role of omitted Fourier modes using the displayed far-field sum.
%
% (b) FINUFFT accuracy.  At nf_full, run at least three values of `nufft_tol`,
%     for example 1e-4, 1e-8, and 1e-12.  Report the same errors and explain
%     any error plateau using the represented Fourier modes, the direct
%     reference, and floating-point arithmetic.  This comparison does not
%     measure `split_tol`, because both calculations use the same split S.
%     Do not infer a complexity or throughput claim from these small runs.
%
% Use `ewald_direct` only as a small-N reference.  Its direct Fourier loop is
% not the fast algorithm being implemented in this worksheet.

%% 6. Submit
%
% Submit one MATLAB script or live script containing:
%
%   1. The local direct routine and the FINUFFT calculation, with the signs,
%      L^(-3) scaling, zero-mode removal, and self correction visible.
%   2. The two requested kernel/error plots and a table for the two studies.
%   3. The MATLAB, FINUFFT, Chebfun, and ESP companion-code versions used.
%   4. A short answer to each question below.
%
% Questions:
%
%   * Why does type 1 create the inner sum in the far-field formula, while
%     type 2 evaluates the outer sum at particle locations?
%   * Which term would be incorrect if the q_i/(S.C0*rc) correction were
%     omitted, and why?
%   * Which parts of this prototype differ from a complete ESP implementation?
%
% For the last question, mention only the particle-grid spreading/interpolation
% and grid-based Fourier evaluation described in the Day 3 notes.  Do not make
% performance claims that are not measured in your calculation.
