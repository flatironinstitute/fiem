r = 1.3; % source box size
R = 5*r; % closest target evaluation

om = 100.1; % wavenumber
nmax = 1000;
js = besselj(0:nmax,sqrt(2)*r*om);
hs = besselh(0:nmax,R*om);

figure(1)
clf
semilogy(abs(hs.*js), 'k.'); hold on;
semilogy(ones(2,1)*sqrt(2)*r*om + log(1/eps), [1e-16,1], 'r--', 'LineWidth',2);
ylim([1e-16,1])
