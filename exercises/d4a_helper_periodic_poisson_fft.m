function [uper,uhat,f0] = d4a_helper_periodic_poisson_fft(f,L)
%D4A_HELPER_PERIODIC_POISSON_FFT Periodic FFT particular solution on [-L,L)^2.
%
%   [UPER,UHAT,F0] = D4A_HELPER_PERIODIC_POISSON_FFT(F,L) solves the mean-zero
%   portion of Delta u = f on the N-by-N uniform grid [-L,L)^2.  F must be
%   stored in the same MESHGRID ordering used by the Day 4 exercise.  UHAT is
%   in FINUFFT's default centered mode order and F0 is the discrete mean.
%
%   The complete particular solution is
%
%       u_p(x1,x2) = u_per(x1,x2) + f0*x1^2/2.
%
%   The second term has to be added separately at boundary targets because it
%   is not periodic.  This small course rewrite follows the centered-mode FFT
%   convention in the Day 4 source derivation and has no private-code
%   dependency.

if ndims(f) ~= 2 || size(f,1) ~= size(f,2)
    error('d4a_helper_periodic_poisson_fft:Grid', ...
        'F must be a square two-dimensional array.');
end
if ~isscalar(L) || L <= 0
    error('d4a_helper_periodic_poisson_fft:Box', 'L must be positive.');
end

N = size(f,1);
if mod(N,2) ~= 0
    error('d4a_helper_periodic_poisson_fft:EvenGrid', ...
        'Use an even number of grid points in each direction.');
end

mid = N/2 + 1;
kmax = pi*N/(2*L);
dk = 2*kmax/N;
k = -kmax + (0:N-1)*dk;
[kx,ky] = meshgrid(k,k);
k2 = kx.^2 + ky.^2;

fhat = ifftshift(fft2(fftshift(f)));
f0 = real(fhat(mid,mid))/N^2;
greenhat = zeros(N,N);
nonzero = k2 ~= 0;
greenhat(nonzero) = -1./k2(nonzero);
uhat = fhat.*greenhat;
uper = real(ifftshift(ifft2(fftshift(uhat))));
end
