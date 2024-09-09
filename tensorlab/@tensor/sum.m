function [ t ] = sum ( t , varargin )
% 
% SUM  Sum of elements.
%    S = SUM(X) is the sum of the elements of the matrix X. S is a row vector with the sum over each column. For N-D arrays, SUM(X) operates along the first
%    non-singleton dimension.
% 
%    S = SUM(X,DIM) sums along the dimension DIM. 
% 
%    Examples:
%    Let  e = cartesianbasis2d('ex','ey')
%
%    If   A = [ e(1) e(2)
%               e(1) e(2)
%               e(2) e(1) ];
% 
%    then sum(X,1) is [ 2.0000*ex + 1.0000*ey, 1.0000*ex + 2.0000*ey ]
%
%    and  sum(X,2) is [ 1.0000*ex + 1.0000*ey
%                       1.0000*ex + 1.0000*ey
%                       1.0000*ex + 1.0000*ey ]
%

% set defaults
dim = 1;

% change dim if t is a row vector
if     size(t,1)==1, dim=2;
elseif size(t,2)==1, dim=1;
end

% check/read the input
nar=size(varargin,2);
if     nar==1, dim=varargin{1};
elseif nar~=0, error('Too many input arguments');
end

% check the dimension
if     ~isfloat(dim),          error('DIM must be a number');
elseif ~all(size(dim)==[1 1]), error('DIM must be a single number'); 
end

% reconstruct input
b  = t.basis;                        % basis
sm = t.size;                         % size of the matrix
st = size(b,1)*ones(1, t.order);     % counter accounting for dimension and order
c  = reshape(t.components, [sm st]); % reshape matrix

% do the summation
c  = sum(c,dim);       % summation
ns = [1 1];            % new matrix size
nd = setdiff(1:2,dim); % other direction than for the summation
ns(nd) = sm(nd);       % copy the relevant matrix size

% reconstruct output
t.size       = ns;     % size
t.components = c(:)';  % components

end
