function [ X ] = zeros ( order, e , varargin )

% ZEROS  Zeros array of tensors
%   X = ZEROS(order,basis,m,n)
%    defines a m x n matrix of zero tensor objects. The tensor objects are
%    defined by the order and the basis given. I.e. if order is 1, we get a
%    m x n matrix of vectors.
% 
%   X = ZEROS(order,basis,[m n p ...])
%    defines a m x n x p x ... multidimensional array of tensor objects
%    defined by the order and the basis
% 
%    For example
%    e = cartesianbasis2d;
%    K = zeros(2,e,10,10);
%    yields a 10 x 10 matrix of 2D zero tensors, with basis e.
% 
%    Consult <a href = "matlab:help zeros">help zeros</a> for more information on the general Matlab zeros

if ~isfloat(order)
    error('Tensor:input','The order must be a number')
    
elseif ~isa( e ,'tensor')
    error('Tensor:input','The basis must be a tensor')
    
elseif ~checkbasis( e )
    error('Tensor:basis','Incorrect basis input')
    
else
    for i = 1:size(varargin,2)
        if ~isfloat(varargin{i})
            error('Tensor:input','Matrix dimensions must be numbers')
        end
    end
    
end

if nargin == 2
    mdim = [1 1];
else
    mdim = [varargin{:}];
end

dim  = size(e.basis,1);

c = zeros( [ mdim dim*ones(1,order) ] );

X.basis = e.basis;
X.order = order;
X.size  = mdim;
X.components = c(:);

X = tensor(X);

if X.order == 0
    X = tensor2double(X);
end

end