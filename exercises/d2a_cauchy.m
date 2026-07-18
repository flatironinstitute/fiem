%% Define family of functions
clear;
close('all')

% This code is meant to demonstrate the failure mode
% of generalized chebyshev procedures. The cauchy kernel
% is nearly hyper singular and thus causes serious difficulties
% and it is difficult to obtain an efficient compressed representation
% for a family of functions which is better than the naive 
% representation
x = chebfun('x');



x0 = -1:0.05:1;
y0 = 0.1:0.025:0.5;

[xx,yy] = meshgrid(x0,y0);
xx = xx(:);
yy = yy(:);
A = [];
tic
for i=1:length(xx)
   A = [A chebfun(1/((x-xx(i))^2 + yy(i)^2))];
end
toc;

[~,nfun] = size(A);

rng(1);
%% generate a random matrix 
B = randn(nfun);

%% Use svd instead of QR, slight abuse of notation to reuse older codes
% Orthogonalize family of functions
[Q,R,V] = svd(A);


%% Use diagonal entries of R as an estimate for rank of family of functions
% Plot diagonal entries of R
figure(1)
clf
semilogy(abs(diag(R)),'k.');

%% Plot the first 6 orthogonal polynomials
figure(2)
clf
plot(Q(:,1:6),'LineWidth',1.6);



%% find number of basis functions which would approximate the full
% family of functions to given precision;

iind_min = min(find(abs(diag(R))<1e-10*abs(R(1,1)))); % don't try and use more,
                                           % than sqrt(eps_mach). Can run
                                           % into stability issues.
                                           % Best to do this in extended 
                                           % precision as a one time calc


% Find good interpolation nodes
iind = min(iind_min+1,nfun);

figure(4)
clf
plot(Q(:,iind+1),'LineWidth',1.6);
ff = Q(:,iind+1); hold on;
rr = roots(ff);
plot(rr,ff(rr),'k.');


ffuns = Q(:,1:length(rr));
B = ffuns(rr);

%% Test accuracy on random function in family 

x0 = 0.01;
y0 = 0.11;

fex = chebfun(y0^2./((x-x0).^2 + y0.^2));

% Construct function samples
fvals = fex(rr);

% Apply values -> coefs matrix, note that inv(B) can be precomputed
% and stored since it is small
coefs = B\fvals;
finterp = ffuns*coefs;

figure(5)
plot(fex,'r-'); hold on; plot(finterp,'k-')

err_interp = norm(fex-finterp);
fprintf('error in interpolating an off-basis function=%d\n',err_interp);






