function t = product(ntc, varargin)

% PRODUCT Product of (matrices of) tensor objects
%   PRODUCT(n,A,B) computes the n-fold inner product of tensors A and B. A
%   and B can also be vectors or matrices of scalars, vectors or tensors.
%   For n=0 the dyadic product results, for n=1 the single inner product,
%   for n=2 the double inner product, etc.
%
%   PRODUCT(n,A,B,C,...) computes the n-fold inner product of an arbitrary
%   number of (matrices of) tensors, where the factors A, B, C, etc. need
%   not be of the same tensorial order.
%
%   PRODUCT is normally called through one of the the methods DYADIC, DOT,
%   DDOT, DDDOT, etc.

% check number and type of input arguments
if nargin < 3
    error('Not enough input arguments.')
end
if ~isnumeric(ntc) || (round(ntc) ~= ntc) || (ntc < 0)
    error('Number of contractions should be nonnegative integer.')
end

% check each tensor in product and determine common basis
b = [];
for i = 1:nargin-1
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
end

% determine sizes of first factor in product and cast it into
% multidimensional array; this array will be used to store all intermediate
% results
t = varargin{1};
if ~isa(t, 'tensor')
    t = double2tensor(t, b);
end
sm = t.size;
st = size(b, 1) * ones(1, t.order);
nm = length(sm);
nt = length(st);
c = reshape(t.components, [sm st]);

% process products one by one
for i = 2:nargin-1
    
    % determine sizes of second factor in product and cast it into
    % multidimensional array
    ti = varargin{i};
    if ~isa(ti,'tensor')
        ti = double2tensor(ti, b);
    end
    smi = ti.size;
    sti = size(b, 1) * ones(1, ti.order);
    nmi = length(smi);
    nti = length(sti);
    ci = reshape(ti.components, [smi sti]);
    
    % check if sufficient number of tensor dimensions left for contraction
    if (nt < ntc) || (nti < ntc)
        error('Number of contractions exceeds tensor order.')
    end

    % if one of the arguments is a 1x1 matrix, no matrix contraction;
    % otherwise check inner matrix dimensions
    if all(sm == 1) || all(smi == 1)
        nmc = 0;
    elseif sm(end) == smi(1)
        nmc = 1;
    else
        error('Inner matrix dimensions must agree.')
    end
    
    % cast multidimensional array of first factor into two-dimensional
    % matrix so that row dimension can be preserved and column dimension is
    % contracted
    c = permute(c, [1:nm-nmc nm+(1:nt-ntc) nm+1-(1:nmc) nm+nt+1-(1:ntc)]);
    c = reshape(c, prod([sm(1:end-nmc) st(1:nt-ntc)]),  prod([sm(nm+1-(1:nmc)) st(nt+1-(1:ntc))]));

    % same for second factor, however here rows will be contracted and
    % columns are preserved
    ci = permute(ci, [1:nmc nmi+(1:ntc) nmc+1:nmi nmi+(ntc+1:nti)]);
    ci = reshape(ci, prod([smi(1:nmc) sti(1:ntc)]), prod([smi(nmc+1:end) sti(ntc+1:end)]));

    % compute actual product
    c = c * ci;
    
    % recast result into multidimensional array with matrix dimensions
    % leading and tensor dimensions following
    c = reshape(c, [sm(1:end-nmc) st(1:end-ntc) smi(nmc+1:end) sti(ntc+1:end)]);
    c = permute(c, [1:nm-nmc nm-nmc+nt-ntc+(1:nmi-nmc) nm-nmc+(1:nt-ntc) nm-nmc+nt-ntc+nmi-nmc+(1:nti-ntc)]);

    % determine new tensor and matrix dimensions and remove superfluous
    % dimensions due to scalar product
    st = [st(1:end-ntc) sti(ntc+1:end)];
    if nmc == 1
        sm = [sm(1:end-nmc) smi(nmc+1:end)];
    else
        if all(sm == 1)
            sm = smi;
        end
        c = reshape(c, [sm st]);
    end
    nm = length(sm);
    nt = length(st);

end

% store result in output tensor
t.order = nt;
t.size = sm;
t.components = c(:);

% turn result into ordinary matrix if tensor character has been lost
if nt == 0
    t = tensor2double(t);
end


end