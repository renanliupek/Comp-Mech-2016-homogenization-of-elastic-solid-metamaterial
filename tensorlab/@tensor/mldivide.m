function [ u ] = mldivide ( varargin )

% MLDIVIDE Backslash or left matrix divide.
% Solves the system K*u = f for u. The order, basis and size of u are
% defined by the order, basis and size of f. The number of contractions is
% determined by the order of K and the order of f.
% 
% Usage: u = K\f
% 
% for the solution method for systems of equation containing matrices of 
% tensor of different order:
% 
% See also CELL/MLDIVIDE

% MLDIVIDE can be used with multiple input arguments to solve a system of
% equations, of which the submatrices are of different order. In this case
% the function should be called by
% 
% mldivide( K(1,1) , K(1,2) , ... , K(1,n) , K(2,1) , ... , K(2,n) , ... ,
% K(n,1) , ... , K(n,n) , f(1,1) , f(2,1) , ... , f(n,1) )
% 
% Here the indices correspond to the different submatrices. Each of this
% matrices can be if size: mi x ni (and does not have to be square by
% definition. As long as the "blown up" matrix is square.

% Determine the dimensions of K and f
n = sqrt(nargin+1/4) - 1/2;

if round(n) ~= n
    error('The number of input arguments does not correspond to K = n x n and f = n x 1')
end

% Dimensions of K and f
Ks = [n n];
Fs = [n 1];

if prod(Ks)+prod(Fs) ~= nargin
    error('The number of input arguments in incorrect')
end

% Positions of K and f in the input
K  = 1:prod(Ks);
F  = prod(Ks)+(1:prod(Fs));

% Check each tensor and determine common basis
b = [];
for i = [K F]
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

if Ks(1) ~= Fs(1)
    error('dimensions do not match');
end

% Define empty cell to store the "blown up" matrices in
KC = {};
FC = {};

% Loop over all indices of K and f
for ii = 1:Ks(1)
 for jj = 1:Ks(2)
  % determine the submatrices of interest. At each time we must no both k
  % and f to determine the number of contractions
  k = varargin{ K( (ii-1)*Ks(2) + jj ) };
  f = varargin{ F( (ii-1)*Fs(2) + 1  ) };
          
  % make sure that k is always a tensor and determine constants
  if ~isa(k, 'tensor')
    k = double2tensor(k, b);
  end
  smk = k.size;                         % matrix dimensions
  stk = size(b, 1) * ones(1, k.order);  % order and dimension (1D/2D/3D)
  nmk = length(smk);                    % order of materix
  ntk = length(stk);                    % order of tensor
  
  % cast k into a multi-dimensional array
  ck = reshape(k.components, [smk stk]);
   
  % make sure that f is always a tensor and determine constants
  if ~isa(f,'tensor')
   f = double2tensor(f, b);
  end
  smf = f.size;
  stf = size(b, 1) * ones(1, f.order);
  nmf = length(smf);
  ntf = length(stf);
  
  % cast f into a multi-dimensional array
  cf = reshape(f.components, [smf stf]);
    
  % determine the number of constraction
  if k.order > f.order
   ntc = k.order - f.order;
  elseif k.order == f.order
   ntc = 0;
  else
   error('K cannot be of lower order than f')
  end
    
  % check if k is of sufficient order to employ the number of contractions
  if ntk < ntc
   error('Number of contractions exceeds tensor order.')
  end

  % if one of the arguments is a 1x1 matrix, no matrix contraction;
  % otherwise check inner matrix dimensions
  if all(smk == 1) || all(smf == 1)
   nmc = 0;
  elseif smk(end) == smf(1)
   nmc = 1;
%   else
%    error('Inner matrix dimensions must agree.')
  end
    
  if k.order > 0
   % cast multidimensional array of first factor into two-dimensional
   % matrix so that row dimension can be preserved and column dimension is
   % contracted
   ck = permute(ck, [1:nmk-nmc nmk+(1:ntk-ntc) nmk+1-(1:nmc) nmk+ntk+1-(1:ntc)]);
   ck = reshape(ck, prod([smk(1:end-nmc) stk(1:ntk-ntc)]),  prod([smk(nmk+1-(1:nmc)) stk(ntk+1-(1:ntc))]));
  end
    
  % add the submatrix to the cell containing all submatrices
  KC{ii,jj} = ck;

  if jj == 1
   if f.order > 0
    % same for second factor, however here rows will be contracted and
    % columns are preserved
    cf = permute(cf, [1:nmc nmf+(1:ntc) nmc+1:nmf nmf+(ntc+1:ntf)]);
    cf = reshape(cf, prod([smf(1:nmc) stf(1:ntc)]), prod([smf(nmc+1:end) stf(ntc+1:end)]));
   end
   
   % add the submatrix to the cell containing all submatrices
   FC{ii,1} = cf;
        
  end
 end
end

% Find the dimensions of the "blown up" matrix. This is an additive
% process, and therefore the dimensions have to start at zero.
KMs = [0 0];
FMs = [0 0];

for ii = 1:Ks(1)
 for jj = 1:Ks(2)
  if size(KC{ii,jj},1) ~= size(FC{ii,1},1)
   error('Incorrect dimensions')
  end
  
  if ii == 1
   KMs(2) = KMs(2) + size(KC{ii,jj},2);
   if jj == 1
    FMs(2) = FMs(2) + size(FC{ii,jj},2);
   end
  end
 
 end
 KMs(1) = KMs(1) + size(KC{ii,1},1);
 FMs(1) = FMs(1) + size(FC{ii,1},1);
end

% error is K is not square
if KMs(1) ~= KMs(2)
    error('The order of K is not consistent')
end

% Empty "blown up" matrix
KM = zeros(KMs);
FM = zeros(FMs);

% Assembly of "blown up" matrix
i=0;j=0;
for ii = 1:Ks(1)
 for jj = 1:Ks(2)
  h = size(KC{ii,jj},1);
  w = size(KC{ii,jj},2);
  KM(i+(1:h),j+(1:w)) = KC{ii,jj};
  j = j+w;
 end
 FM(i+(1:h),:) = FC{ii,1};
 i=i+h;
 j=0;
end

% Solve the "blown up" system of equations using Matlab's standard MLDIVIDE
% function (i.e. matrix division of KM and FM)
UM = KM\FM;

% Rebuild the solution as a system of doubles and tensor objects
U = {};

i=0;
for ii = 1:size(FC,1)
   
 % Get the subentry from f (we have assumed f and u of equal order and 
 % basis)
 u = varargin{ F( (ii-1)*Fs(2) + 1  ) };
 
 % Make sure that u is a tensor object (even for a complete scalar entry)
 % to do all the operations
 if ~isa(u,'tensor')
  u = double2tensor(u, b);
 end
 
 % Get the relevant dimensions of u
 smu = u.size;
 stu = size(b, 1) * ones(1, u.order);
 
 % Get the subentry from the "blown up" matrix UM
 c = UM(i+(1:size(FC{ii,1},1)),size(FC{ii,1},2));
 i = i + size(FC{ii,1},1);
 
 % Cast c into a materix suitable to construct the tensor object
 c = reshape(c,[prod([stu smu]) 1]);

 % Restore u as a tensor
 u.basis = b;
 u.size = smu;
 u.components = c;

 % turn result into ordinary matrix if tensor character has been lost
 if u.order == 0
  u = tensor2double(u);
 end
 
 % Assembly into the the system of results
 U{ii,1} = u;
 
end

% If U is a single entry (column of scalars/tensors) give the result as a
% matrix object. Else return the result as a cell.
if size(U,1) == 1 && size(U,2) == 1
    u = U{:};
else
    u = U;
end

% End of MLDIVIDE
end

