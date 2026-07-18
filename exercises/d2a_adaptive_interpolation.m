%WORKSHEET2_adaptive_function_interpolation
%
%
% Requirements:
%   - This worksheet will use the chunkIE package to demonstrate the ideas.
%   Please download and install the package according to the directions
%   here: https://chunkie.readthedocs.io/en/latest/getchunkie.html
%   - chebfun: For high accuracy representation of the function
%
% Scope:
%   - Adaptive discretizations for nearly singular functions
%
clear
close('all')
% This code steps through how an adaptive Legendre interpolator would
% proceed

% try out different values of p, and different functions f below




p = 16;
eta = 0.01;

f = @(x) eta.^2./(x.^2 + eta.^2);
fex = chebfun(f);


% This is the main part of the code
nlmax = 20; % Maximum level of refinement
[x,w,u,~] = lege.exps(p);
x = (x+1)/2;
xequi = 0:1/p:(1-1/p);

nbmax = 1000;
xnodes = zeros(p,nbmax);
iboxuse = zeros(1,nbmax);
iproclist = zeros(1,nbmax);
fvals = zeros(p,nbmax);

fcoeffs = zeros(p,nbmax);



endpts = zeros(2,nbmax);
endpts(1,1) = -1;
endpts(2,1) = 1;

a = endpts(1,1);
b = endpts(2,1);
xnodes(:,1) = a + (b-a)*x;
fvals(:,1) = f(xnodes(:,1));
fcoeffs(:,1) = u*fvals(:,1); 
iboxuse(1) = 1;
iproclist(1) = 1;
ntot = 1;

tol = 1e-12;
for ilev = 1:nlmax
    jboxes = find(iproclist == 1);
    if(length(jboxes) == 0), break; end
    
    iproclist(jboxes) = 0;
    for jbox = jboxes
        
        
        a = endpts(1,jbox); b = endpts(2,jbox);
        err_test = (abs(fcoeffs(end,jbox)) + abs(fcoeffs(end-1,jbox)))*(b-a)/2;
        
        if(err_test > tol) 
            iboxuse(jbox) = 0;
            box_mid = mean(endpts(1:2,jbox));
        
            a = endpts(1,jbox);
            b = box_mid;
        
            ncur = ntot+1;
            endpts(1,ncur) = a;
            endpts(2,ncur) = b;
            xnodes(:,ncur) = a + (b-a)*x;
            fvals(:,ncur) = f(xnodes(:,ncur));
            fcoeffs(:,ncur) = u*fvals(:,ncur);
            iboxuse(ncur) = 1;
            iproclist(ncur) = 1;

            ncur = ntot+2;
            a = box_mid;
            b = endpts(2,jbox);
        
            endpts(1,ncur) = a;
            endpts(2,ncur) = b;
            xnodes(:,ncur) = a + (b-a)*x;
            fvals(:,ncur) = f(xnodes(:,ncur));
            fcoeffs(:,ncur) = u*fvals(:,ncur);
            iboxuse(ncur) = 1;
            iproclist(ncur) = 1;
        
            ntot = ntot+2;
        end
        
    end
    xplot = xnodes(:,iboxuse==1);
    fplot = fvals(:,iboxuse==1);
    fcoeffsplot = fcoeffs(:,iboxuse==1);


    figure(1)
    clf
    subplot(2,1,1);
    plot(fex,'k-'); hold on;
    for i=1:size(xplot,2)
        plot(xplot, fplot,'r-'); 
    end
    xlim([-1.1,1.1]);
    ylim([-0.2,1.2]);
    abuni = [unique(endpts(1,iboxuse==1)) 1];
    plot([-1,1],[-0.15 -0.15],'k-','LineWidth',1.6);
    plot(abuni,-0.15*ones(length(abuni)),'r.','MarkerSize',20)
    axis equal;
    xl = xlim;
    
    subplot(2,1,2)
    abuse = endpts(1:2,iboxuse==1);
    for i=1:size(xplot,2)
        a = abuse(1,i); b = abuse(2,i);
        xuse = a + xequi*(b-a);
        semilogy(xuse,fcoeffsplot(:,i),'k.'); hold on;
    end
    ylim([1e-18,10]);
    xlim(xl);
    pause    
end

fprintf('Number of points in global chebyshev exp=%d\n',...
    length(fex));
nadap=length(xplot(:));
fprintf('Number of points in adaptive Legendre exp=%d\n',...
    nadap);




%%  Exercise 1:
%  Using your spectral Legendre interpolator in (Exercise 2) of the
%  function_interpolator worksheet, compare how well does the adaptive
%  interpolator do as compared to the global one as a function of k
%
%  when eps = 1e-12, for
%  a) log(|x - (1 + 1/k)|)
%  b) sin(pi*k*x)
%


%% Exercise 2:
%  Choosing p based on tolerance. The choice of number of points you
%  use per panel in the adaptive interpolator can be optimized to 
%  minimize the total number of points as a function of tolerance. 
%  A usual rule of thumb in the community is p = \log_{10} (1/\eps) + 1
%
%  Consider the Runge function given in this worksheet. For tolerances
%  of 1e-6, 1e-8, 1e-10, and 1e-12, plot the total number of points (nadap) as
%  a function of p. Does the above trend hold true? 
%  
