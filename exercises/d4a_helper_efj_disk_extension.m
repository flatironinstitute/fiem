function [fext,info] = d4a_helper_efj_disk_extension(X,Y,source,norder,a,r0,r1)
%D4A_HELPER_EFJ_DISK_EXTENSION EFJ extension of a scalar function from the unit disk.
%
%   [FEXT,INFO] = D4A_HELPER_EFJ_DISK_EXTENSION(X,Y,SOURCE,NORDER,A,R0,R1)
%   extends SOURCE from the unit disk to the Cartesian grid specified by X,Y.
%   X and Y must be same-sized arrays, such as those returned by MESHGRID.
%   SOURCE accepts two same-sized arrays x,y and returns its values f(x,y).
%
%   The extension is the finite EFJ formula from the Day 4 morning notes:
%
%     E[f](y + d*n_y) = Phi(d) sum_j w_j f(y - t_j*d*n_y).
%
%   On the unit circle, y=n_y=(x,y)/sqrt(x^2+y^2).  The nodes t_j are the
%   scaled Chebyshev nodes and the weights are the explicit Lagrange-basis
%   weights.  Phi is a C-infinity cutoff that equals one for d <= R0 and zero
%   for d >= R1.  This helper is a small course rewrite of those formulas; it
%   does not depend on the private function-extension repository.

if ~isequal(size(X),size(Y))
    error('d4a_helper_efj_disk_extension:GridSize', ...
        'X and Y must have the same size.');
end
if ~isa(source,'function_handle')
    error('d4a_helper_efj_disk_extension:Source', ...
        'SOURCE must be a function handle.');
end
if ~isscalar(norder) || norder < 0 || norder ~= floor(norder)
    error('d4a_helper_efj_disk_extension:Order', ...
        'NORDER must be a nonnegative integer.');
end
if ~isscalar(a) || a <= 0
    error('d4a_helper_efj_disk_extension:NodeRange', ...
        'A must be positive.');
end
if ~isscalar(r0) || ~isscalar(r1) || r0 < 0 || r1 <= r0
    error('d4a_helper_efj_disk_extension:Window', ...
        'Require 0 <= R0 < R1.');
end
if a*r1 > 1 + 64*eps
    error('d4a_helper_efj_disk_extension:SampleOutsideDisk', ...
        'Require A*R1 <= 1 so every EFJ sample remains in the unit disk.');
end

[t,w] = efj_nodes_weights(norder,a);
r = hypot(X,Y);
inside = r <= 1;
d = r - 1;
extension = (d > 0) & (d < r1);

fext = zeros(size(X));
fext(inside) = source(X(inside),Y(inside));

if any(extension(:))
    dext = d(extension);
    nx = X(extension)./r(extension);
    ny = Y(extension)./r(extension);
    value = zeros(size(dext));

    for j = 1:numel(t)
        sample_radius = 1 - t(j)*dext;
        value = value + w(j)*source(sample_radius.*nx,sample_radius.*ny);
    end

    fext(extension) = smooth_cutoff(dext,r0,r1).*value;
end

info = struct();
info.inside = inside;
info.extension = extension;
info.distance = d;
info.nodes = t;
info.weights = w;
info.minimum_sample_radius = 1 - max(t)*r1;
end

function [t,w] = efj_nodes_weights(norder,a)
if norder == 0
    t = 0;
    w = 1;
    return
end

j = (0:norder).';
t = 0.5*a*(1-cos(j*pi/norder));
x0 = 1 + 2/a;
s = sqrt(x0^2 - 1);
tnp = (1+a)*((x0+s)^norder-(x0-s)^norder)/s;
w = (-1).^j*tnp./(norder*a*(1+t));
w(1) = w(1)/2;
w(end) = w(end)/2;
end

function value = smooth_cutoff(d,r0,r1)
value = zeros(size(d));
value(d <= r0) = 1;
transition = (d > r0) & (d < r1);
if any(transition(:))
    s = (d(transition)-r0)/(r1-r0);
    left = exp(-1./s);
    right = exp(-1./(1-s));
    value(transition) = right./(left+right);
end
end
