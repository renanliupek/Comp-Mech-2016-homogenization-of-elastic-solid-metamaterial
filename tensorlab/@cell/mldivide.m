function [ u ] = mldivide ( varargin )

% MLDIVIDE Backslash or left matrix divide.
% Solves the system K*u = f for u. The order, basis and size of u are
% defined by the order, basis and size of f. The number of contractions is
% determined by the order of K and the order of f.
% 
% Usage the system of equation is formed by (e.g. K = 2x2)
%  K = { K(1,1) , K(1,2) ; K(2,1) , K(2,2) }
%  f = { f(1,1) ; f(2,1) }
% 
%  u = K\f
% 
% Herin, each component can be a matrix of tensors of any order, or an
% ordinary Matlab matrix of doubles.
% 
% For the solution of a system with a single matrix K of equal order, with 
% a column f of equals order as the right-hand side:
% 
% See also TENSOR/MLDIVIDE 
% 

% Reformulate the input to conform to tensor/mldivide
Ko = varargin{1};
Fo = varargin{2};

K = {};
F = {};
for ii = 1:size(Ko,1)
    K = [K Ko(ii,:)];
    F = [F Fo(ii,:)];
end

u = mldivide(K{:},F{:});

end