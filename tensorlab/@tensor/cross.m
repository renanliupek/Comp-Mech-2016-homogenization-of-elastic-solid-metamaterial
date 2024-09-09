function t = cross(varargin)

% CROSS  Vector/tensor cross product.
%   CROSS(A,B), where A and B are vectors or tensors, returns the cross
%   product AxB of A and B. A and B must have a 3D cartesian basis.
%
%   CROSS(A,B), where A and B are matrices of 3D vectors or tensors,
%   returns the cross product of these matrices. The usual contraction
%   applies to the matrix dimensions of A and B and their inner dimensions
%   must correspond. If A or B is a 1x1 matrix, it is multiplied with all
%   components of the remaining matrix.
%
%   CROSS(A,B,C,...) returns the cross products of tensor-matrices A, B, C,
%   etc.
%
%   See also DYADIC, DOT, DDOT.

if nargin < 2
    error('Not enough input arguments.')
end

% set basis
b = [];
for i = 1:nargin
    ti = varargin{i};
    if isa(ti, 'tensor')
        if isempty(b)
            b = ti.basis;
        elseif size(ti.basis, 1) ~= size(b, 1) || any(any(ti.basis ~= b))
            error('Tensor bases must agree.')
        end
    elseif ~isnumeric(ti)
        error('Multiplication of tensor and %s is undefined.', class(ti))
    end
    if ndims(ti) > 2
        error('Input arguments must be two-dimensional matrices.')
    end
end
if isempty(b)
    error('No tensor object found in product.')
elseif size(b, 1) ~= 3
    error('Tensor basis must be 3D.')
end

% prepare base tensor
t = varargin{1};

% construct transformation tensor
eta.basis = b;
eta.order = 3;
eta.size = [1 1];
eta.components = [0 0 0 0 0 1 0 -1 0 0 0 -1 0 0 0 1 0 0 0 1 0 -1 0 0 0 0 0]';
eta = tensor(eta);

% process individual cross products by rewiting them to inner products
for i = 2:nargin
    ti = varargin{i};
    t = dot(t, eta, ti);
end

end