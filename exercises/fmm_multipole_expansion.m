%WORKSHEET_multipole_expansions
%
%
% Scope:
%   - Multipole expansion lengths for Laplace and Helmholtz
%

boxsize = 0.5;
nsrc = 100;
src = (rand(2,nsrc) -0.5)*boxsize;
zsrc = src(1,:) + 1j*src(2,:);

ctarg = [2*boxsize;0]; % presribe distance between centers/center of target box

ntarg = 200;
targ = sep + (rand(2,ntarg) - 0.5)*boxsize;
ztarg = targ(1,:) + 1j*targ(2,:);

charges = rand(nsrc,1);

pot_ex = 1./(ztarg.' - zsrc)*charges;

% Test form multipole and multipole eval routine for Laplace
n = 50;
alpha = form_lap_multipole(zsrc, charges, n);
pot = eval_lap_multipole(alpha, ztarg, n);

err = norm(pot-pot_ex)./sum(abs(charges(:)));
fprintf('Error in multipole expansion=%d\n',err);


% Test form multipole and multipole eval routine for Helmholtz
zk = 3.5;
n = 50;
rr = abs(ztarg.' - zsrc);
pot_helm_ex = besselh(0,zk*rr)*charges;
alpha_helm = form_helm_multipole(zk, zsrc, charges, n);
pot_helm = eval_helm_multipole(zk, alpha_helm, ztarg, n);

err_helm = norm(pot_helm-pot_helm_ex)./sum(abs(charges(:)));
fprintf('Error in multipole expansion=%d\n',err_helm);



%% Exercise 1:
% Determine the number of terms required to achieve a tolerance \eps. Plot
% N as a function of \eps, verify that N agrees with your estimate of it
% based on the separation chosen. b) Is N a function of the boxsize
% 
% Repeat the above exercise for the Helmholtz case. Does the dependence
% on boxsize agree with the expectation



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


function alpha = form_helm_multipole(zk, zsrc, charges, n)
    alpha = zeros(2*n+1,1);
    thet = atan2(imag(zsrc), real(zsrc));
    for l=-n:n
        js = besselj(l,zk*abs(zsrc(:))).*exp(-1j*l*thet(:));
        alpha(l+1+n) = sum(js.*charges);
    end
end

function pot = eval_helm_multipole(zk, alpha, ztarg, n)
    nt = length(ztarg);
    pot = zeros(nt,1);
    thet = atan2(imag(ztarg), real(ztarg));
    for l=-n:n
        hs = besselh(l,zk*abs(ztarg(:))).*exp(1j*l*thet(:));
        pot = pot + alpha(l+1+n)*hs(:);
    end
end
