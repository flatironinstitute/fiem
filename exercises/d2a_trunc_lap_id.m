%% Define family of functions
clear
close('all')


npts = 3000;
[x,w] = chebpts(npts);
x = (x + 1)/2;
w = w/2;


lj = 0.01:0.001:1;
A = exp(-x./lj).*sqrt(w.');
[~,nfun] = size(A);

tol = 1e-7;

[sk, ~, ~] = id(A', tol);

rng(1);

[Qfun, Rfun] = qr(A, 'econ');

%% Find row indices using pivoted qr
%% Use diagonal entries of R as an estimate for rank of family of functions
iind_min = min(find(diag(abs(Rfun))<tol*abs(Rfun(1,1)))); % don't try and use more,
                                           % than sqrt(eps_mach). Can run
                                           % into stability issues.
                                           % Best to do this in extended 
                                           % precision as a one time calc
iind = min(iind_min+1,nfun);

[~,~,P] = qr(Qfun(:,1:iind)', 'econ');

[p, ~] = find(P');
sk2 = p(1:iind);
%%
figure(1)
clf
semilogy(diag(abs(Rfun)),'k.');


%% Plot the first 6 orthogonal functions
Q = Qfun./sqrt(w.');
figure(2)
clf
plot(x,Q(:,1:6),'LineWidth',1.6);

% Plot the same 6 orthogonal polynomials but with a log scale in x
% to see the structure close to the origin
figure(3)
clf
semilogx(x,Q(:,1:6),'LineWidth',1.6);


%% find number of basis functions which would approximate the full
% family of functions to given precision;


rr = x(sk);
rr2 = x(sk2);


B = Qfun(sk,1:length(sk));
B2 = Qfun(sk2,1:iind);

Bmat = Q(sk2,1:iind);

%% Test accuracy on random function in family 


lj_test = 0.01145;
fex = @(x) exp(-x/lj_test);

% Construct function samples
fvals = fex(rr2);

% Apply values -> coefs matrix, note that inv(B) can be precomputed
% and stored since it is small
coefs = Bmat\fvals;
finterp = Q(:,1:iind)*coefs;

figure(5)
plot(x,fex(x),'r-'); hold on; plot(x,finterp,'k-')

err_interp = norm(fex(x)-finterp);
fprintf('error in interpolating an off-basis function=%d\n',err_interp);


%% (Exercise 1) 
% Repeat the same set of exercises but for the functions x^{\mu} for \mu \in [0.5, 40]
% Note that you may need to use exp(mu \log(x)) using chebfun to get the test to work
% What is the numerical rank of this family of functions? 
%

%% (Exercise 2)
% Repeat the above exercise but for \mu \in [-0.5, 40]. What is the rank in this case?  




