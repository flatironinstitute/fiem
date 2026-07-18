%WORKSHEET1_function_interpolation
%
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%
% Scope:
%   - Adaptive discretizations for nearly singular functions
%

%% (Section 1) Legendre nodes 
%
% Legendre nodes are points in [-1,1] which are good for interpolating and
% integrating functions defined on [-1,1]
% 

%% 
% Legendre nodes and polynomials

close all; clearvars; clc;

% get some quantities for working with Legendre polynomials of order k

k = 16;
[x,w,u,v] = lege.exps(k);

% x -> Legendre nodes
% w -> Legendre weights
% u -> kxk matrix maps values at points x to Legendre series coefficients
%      of polynomial interpolant
% v -> kxk matrix evaluates Legendre series at points x. note v = inv(u)

figure(1)
clf
plot(x,zeros(size(x)),'rx')
hold on 

% can get Legendre polynomial values at various points from lege.pols

xviz = linspace(-1,1,200);
pviz = lege.pols(xviz(:),k);
plot(xviz,pviz(2,:))
plot(xviz,pviz(k/2,:))
plot(xviz,pviz(k+1,:))
legend('','P_1','P_{k/2}','P_{k}')

%%
% we can find the interpolant of f = sin(x) and look at the coefficients


% plot sin(n(x-0.1)^2) at 200 equispaced points in [-1,1] and plot values at
% Legendre nodes

n = 5;
ffun = @(t) sin(n*(t-0.1).^2);

fvals = ffun(x);
fviz = ffun(xviz);

figure(2)
clf
plot(xviz,fviz,'k-')
hold on
plot(x,fvals,'bo')

% compute values of interpolant at these points and compare to interpolant

fcoefs = u*fvals;
interpmat = pviz(1:k,:).';
finterp = interpmat*fcoefs;
plot(xviz,finterp,'b--')

% plot coefficients on a semilogy plot

figure(3)
clf
semilogy(abs(fcoefs),'b-o')

%%
% the integral obtained from sum(fvals.*w) is accurate for small n

fint_true = integral(ffun,-1,1,"RelTol",1e-15);
fint_approx = sum(fvals.*w);
abs(fint_true-fint_approx)

%% (Exercise 1)
%
% - use the Legendre nodes and weights to compute the integral of a 
% the function f = log(|t-1.1|) over [-1,1]
% - compare to an accurate integral obtained from MATLAB's integral
% function
% - see how the error depends on the order k 

%% (Exercise 2)
%
%  Write a routine that adaptively determines the order of legendre
%  expansion required to resolve the function f. You can use two different
%  monitor functions to determine if the function is resolved to a
%  tolerance of \eps
%
%  Option 1: Test the interpolation error at a finer set of nodes
%  Option 2: Use the last two coefficients of the Legendre function
%
%  Using your spectral function discretizer, determine the order of
%  Legendre expansion (q) required for obtaining tolerance \eps as a function
%  of k for the following family of functions
%
%  a) log(|x - (1 + 1/k)|)
%  b) sin(pi*k*x)
%
%  vary k \in [1,1000]
%
%  How does the legendre order q depend on k for a fixed tolerance?

