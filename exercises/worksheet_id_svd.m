%WORKSHEET_id_svd
%
% 
% Scope:
%   - Multipole expansions and comparison to svd
% Related worksheet: fmm_multipole_expansion.m

% Consider the kernel in the Laplace FMM again. Let p denote the number 
% of terms determined in worksheet_multipole_expansion that you derived. 

% Create sources distributed uniformly in the box of diameter boxsize
boxsize = 1.0;
nsrc = 100;
src = (rand(2,nsrc) -0.5)*boxsize;
zsrc = src(1,:) + 1j*src(2,:);

ctarg = [2*boxsize;0]; % presribe distance between centers/center of target box

% Create targets distributed uniformly in box of diameter boxsize
% but centered at ctarg
ntarg = 200;
targ = ctarg + (rand(2,ntarg) - 0.5)*boxsize;
ztarg = targ(1,:) + 1j*targ(2,:);



% Form the kernel matrix
A = 1./(ztarg.' - zsrc);

% Apply the kernel matrix to a given vector
charges = rand(nsrc,1);
pot_ex = A*charges;

% Compute and evaluate the truncated multipole expansion
n = 50; % set the number of terms in the multipole expansion here
alpha = form_lap_multipole(zsrc, charges, n);
pot = eval_lap_multipole(alpha, ztarg, n);

err = norm(pot-pot_ex)./sum(abs(charges(:)));
fprintf('Error in multipole expansion=%d\n',err);

% Compute the svd of the kernel matrix
ss = svd(A);

figure(1)
clf()
semilogy(ss, 'b-o', 'MarkerSize', 4);
xlabel('index'); ylabel('singular value');
title('Singular values of A (well-separated Cauchy kernel)');
grid on;

%% Exercise 1 
% a) How does the best rank using svd compare to the number of terms
% in the multipole expansion?
%
% b) Redo the above example, but in this case, suppose the sources are
% uniformly distributed in the circle of radius boxsize, and that the
% targets are uniformly distributed in the annulus between radii [R,R+1]
%
% How does the svd rank compare to the multipole rank in this case?
%
% c) As with the multipole expansion case, how does the numerical rank
%    depend on the separation




%% Interpolative decomposition

tols = [1e-2, 1e-4, 1e-8, 1e-12];
fprintf('%8s  %6s  %12s\n', 'tol', 'rank k', '||A-A_ID||/||A||');
for tol = tols
    [sk, rd, T] = id(A, tol);      % FLAM's id function
    k = length(sk);
    % Reconstruct: A(:,rd) ≈ A(:,sk)*T, so full approx is K with rd cols replaced
    A_approx = A;
    A_approx(:, rd) = A(:, sk) * T;
    rel_err = norm(A - A_approx) / norm(A);
    fprintf('%8.0e  %6d  %12.4e\n', tol, k, rel_err);
end

%% Exercise 2
%
% Compare the errors obtained in the previous factorization to the best
% possible error \sigma_{k+1}, where \sigma_{j} are the singular values 
% of the matrix. 
%
% How does the ratio err/\sigma_{k+1} grow as you proportionately increase
% the number of sources and targets? How does it compare to the theoretical
% result of err <= (1+ \sqrt{k(N-k}) \sigma_{k+1}?
%

%% Visualizing the interpolative decomposition
[sk, rd, T] = id(A, 1e-8);

figure; hold on;
plot(real(zsrc(rd)), imag(zsrc(rd)), 'b.', 'MarkerSize', 10, 'DisplayName', 'redundant');
plot(real(zsrc(sk)), imag(zsrc(sk)), 'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'skeleton');
plot(real(ztarg), imag(ztarg), 'k.', 'MarkerSize', 6, 'DisplayName', 'targets');
axis equal
legend; 
title(sprintf('ID skeleton (k=%d) vs redundant (N=%d) source points', length(sk), nsrc));
xlabel('position');


%% The power of randomization. 
%
% FLAM's id under the hood uses a rank-revealing QR algorithm to compute
% the skeleton and redundant points (columns). This computation can get
% increasingly expensive if the size of the system increases.
%
%
%% Exercise 3
% Set the number of sources=number of targets = N.
% a) Time the FLAM id routine as a function of N, how does the time required
% grow?
%
% b) Implement the following randomized algorithm, let W be a random 
% N \times (k+10) matrix where k is the known rank (this can be adaptively 
% determined). Compute B = A*W. Compute the id of B. Use the skeleton
% points of B as the skeleton for A, and you may use the same interpolation
% matrix T to extract the remaining columns of A. 
%     i) What is the error in the reconstruction of A using the randomized
%        id as compared to id in FLAM?
%     ii) Time the randomized id routine as a function of N (excluding the
%         random matrix generation time) and see how it compares to 
%         the id of FLAM.

function alpha = form_lap_multipole(zsrc, charges, n)
    alpha = zeros(n,1);
    nsrc = length(zsrc);
    zpow = ones(nsrc,1);
    for l=0:n-1
        alpha(l+1) = sum(zpow.*charges);
        zpow = zpow.*zsrc(:);
    end
end

function pot = eval_lap_multipole(alpha, ztarg, n)
    zpow = 1./ztarg(:);
    zfac = 1./ztarg(:);
    nt = length(ztarg);
    pot = zeros(nt,1);
    for l=0:n-1
        pot = pot + alpha(l+1)*zpow(:);
        zpow = zpow.*zfac;
    end
end

