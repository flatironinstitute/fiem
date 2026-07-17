%D1M_EXTERIOR_LAPLACE_BIE Day 1 morning exercise.
%
% Exterior Dirichlet problem for the two-dimensional Laplace equation
%
% This exercise uses chunkIE.  It is deliberately a short numerical
% investigation rather than a programming project.  You may use AI tools
% to help write MATLAB code, but you are responsible for deriving the two
% equations, running the code, and checking every reported result.
%
% Required software
% -----------------
% Install chunkIE and run its startup script.  Its installation instructions
% are at https://chunkie.readthedocs.io/en/latest/getchunkie.html.
% chunkIE uses FMM2D when the relevant accelerated evaluators are available.
% Also install FINUFFT and FLAM for later exercises, but they are not needed
% for this exercise.
%
% Before beginning, record the MATLAB version and the git revisions (or
% release versions) of the packages that you use.  Verify that the following
% commands find chunkIE functions:
%
%   which chunkerfunc
%   which chunkermat
%   which chunkerkerneval
%
% Problem
% -------
% Let Gamma be one smooth, simple, closed, positively oriented curve, with
% outward normal n.  Let Omega be its bounded interior.  Solve
%
%   Delta u = 0              in R^2 \ closure(Omega),
%   gamma^+ u = f            on Gamma,
%   u(x) = O(1/|x|)          as |x| -> infinity.
%
% Here gamma^+ is the exterior boundary limit.  Use
%
%   G(x,y) = -(1/(2*pi))*log(|x-y|).
%
% Start with chunkIE's starfish curve.  You may replace it with a favorite
% smooth curve, provided that it has one component, has no self
% intersections, and that you verify the selected source and target points
% are in the intended regions.

%% 1. Manufactured solution and boundary data
%
% Select two distinct points z1 and z2 strictly inside Omega and define
%
%   u_ex(x) = G(x,z1) - G(x,z2).
%
% The singularities are excluded from the exterior problem.  The two
% logarithmic far fields cancel, so u_ex(x) = O(1/|x|).  Set
%
%   f = u_ex|_Gamma.
%
% Use `chunkerinterior` to verify that z1 and z2 are inside the curve, and
% choose a fixed set of evaluation targets that are outside the curve.
% Plot the curve, the two interior source points, and the evaluation targets.
%
% Suggested initial geometry and source locations:
%
%   nch = 8;
%   chnkr = chunkerfuncuni(@(t) starfish(t),nch);
%   z1 = [ 0.15;  0.05];
%   z2 = [-0.20; -0.10];
%
% A convenient MATLAB definition of the manufactured field is
%
%   G = @(x,z) -log(sqrt(sum((x-z).^2,1)))/(2*pi);
%   uex = @(x) G(x,z1) - G(x,z2);
%
% Construct f from `uex(chnkr.r)` after you have created the geometry.

%% 2. Derive two boundary integral equations before coding
%
% Define the single- and double-layer potentials by
%
%   S[sigma](x) = int_Gamma G(x,y) sigma(y) ds_y,
%   D[mu](x)    = int_Gamma dG(x,y)/dn(y) mu(y) ds_y.
%
% (a) First-kind formulation.  Use
%
%       u = S[sigma] + c,
%
%     together with the zero-total-charge condition
%
%       int_Gamma sigma ds = 0.
%
%     Derive the bordered boundary system
%
%       [ S   1 ] [sigma] = [ f ],
%       [ w^T 0 ] [  c  ]   [ 0 ],
%
%     where w contains the quadrature weights.  Explain why zero total
%     charge removes the logarithmic far field and why c should converge to
%     zero for this manufactured solution.
%
% (b) Second-kind formulation.  Use
%
%       u = D[mu] + c.
%
%     With the outward normal and exterior trace convention used in the Day
%     1 morning notes,
%
%       gamma^+ D[mu] = (1/2 I + D) mu.
%
%     The density has a constant nullspace.  Fix it with
%
%       int_Gamma mu ds = 0,
%
%     and derive the corresponding bordered system.  Again, explain why the
%     computed c should be close to zero.
%
% Do not change the sign of the jump term to make a computation look better.
% Check the convention against the Day 1 notes and validate it using the
% manufactured solution.

%% 3. Implement the two Nyström discretizations in chunkIE
%
% Assemble the single- and double-layer matrices with, for example,
%
%   kern_s = kernel('lap','s');
%   kern_d = kernel('lap','d');
%   S = chunkermat(chnkr,kern_s);
%   D = chunkermat(chnkr,kern_d);
%
% Use `chnkr.wts(:)` for w and solve the two bordered linear systems from
% Section 2.  The code you submit should make the following quantities
% explicit:
%
%   * the number of nodes and the panel order used by chunkIE;
%   * the vectors of quadrature weights and boundary data;
%   * the two matrices and their jump/constraint terms; and
%   * the constants c and the residuals of the zero-mean constraints.
%
% Evaluate both potentials at the exterior target points with
% `chunkerkerneval`, then add c.  Do not evaluate a singular or nearly
% singular layer potential with a direct smooth-rule call to `kern.eval`.
% The on-surface matrices and near-target evaluations must use chunkIE's
% native singular/corrected quadrature.
%
% Report, for each formulation,
%
%   max_j |u_num(x_j) - u_ex(x_j)|,
%   ||A rho + c*1 - f||_inf,
%   |w^T rho|, and |c|,
%
% where (A,rho) is (S,sigma) or (1/2*I+D,mu), as appropriate.

%% 4. Resolution and conditioning study
%
% Hold the geometry, z1, z2, target points, and the panel order fixed.  Run
% both formulations for at least five increasing values of `nch`, for example
%
%   nch = [4 6 8 12 16].
%
% For every resolution, tabulate the number of nodes, target error, boundary
% residual, constraint residual, |c|, and a condition-number estimate for
% each bordered linear system.  A dense `cond` calculation is sufficient at
% the modest resolutions in this exercise.
%
% Make one convergence plot and one conditioning plot.  State explicitly
% which matrices and which norm/scaling you used for the condition numbers:
% a discrete Euclidean matrix condition number depends on the chosen
% parametrization and scaling.  Compare the observed behavior of the
% first-kind and second-kind formulations without claiming an asymptotic rate
% that your data do not support.

%% 5. Submit
%
% Submit one MATLAB live script, script, or notebook and a short report that
% contains:
%
%   1. The curve and the manufactured solution, including a statement of why
%      it is an admissible exterior solution.
%   2. The derivation of both bordered equations and the exterior jump sign.
%   3. One geometry/source/target plot, a convergence plot, and a conditioning
%      plot.
%   4. A table of the numerical checks from Section 3.
%   5. A brief comparison of the two formulations and any use of AI tools.
%
% Optional extension
% ------------------
% Repeat the target-error test at points approaching Gamma from the exterior.
% Compare chunkIE's corrected evaluation with a deliberately forced smooth
% evaluation, and explain the observed loss of accuracy in terms of the
% quadrature rather than the boundary integral equation itself.
