function t3lt = l3transpose(t,b,o)

%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% L3TRANSPOSE left tranaspose of a third order tensor.
%
%   t = T2VOIGT(t,b,o)
%   
%   input  t  - tensor
%   input  b  - basis
%   input  o  - order of the tensor
%   output t3lt - tensor
%
%--------------------------------------------------------------------------



if o~=3; error("tensor should be of order 3."); end;

% spatial dimension
dim  = size(b, 1);

% case 2D
if dim==2

% % extract basis component
% e1 = b(1);
% e2 = b(2);


% matrix representation
m = zeros(dim,dim,dim);


% for loop to save t_ijk components in a matrix m
for i=1:dim
    for j=1:dim
       for k=1:dim
           % unit vectors ei,ej,ek
           ei=b(i);
           ej=b(j);
           ek=b(k);
           % component i,j,k of tensor t
           m(i,j,k) = dot(ej,dot(ei,t,ek));
       end
    end
end


%initialize left transposed tensor with zero
t3lt = 0*b(1)*b(1)*b(1);

% for loop to assign components of the left transposed tensor
for i=1:dim
    for j=1:dim
       for k=1:dim
           % unit vectors ei,ej,ek
           ei=b(i);
           ej=b(j);
           ek=b(k);
           % left transposed tensor
           t3lt = t3lt+ m(j,i,k)*ei*ej*ek;
       end
    end
end


% % consistency check
% m3lt = zeros(dim,dim,dim);
% % for loop to save t3lt_ijk components in a matrix m3
% for i=1:dim
%     for j=1:dim
%        for k=1:dim
%            % unit vectors ei,ej,ek
%            ei=b(i);
%            ej=b(j);
%            ek=b(k);
%            % component i,j,k of tensor t
%            m3lt(i,j,k) = dot(ej,dot(ei,t3lt,ek));
%        end
%     end
% end

% case 3D
elseif dim==3
    error('l3transpose only supports 2D.');
end


end



