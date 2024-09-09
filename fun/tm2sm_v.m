function S = tm2sm(t,b,o)
%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% TM2SM  Transform a tensor matrix to an equivalent scalar matrix
%   S = TM2SM(t,b,o)
%   
%   input  t - tensor matrix
%   input  b - basis
%   input  o - order of the tensor matrix component
%   output S - equivalent square scalar matrix S 
%
%--------------------------------------------------------------------------

% Comments (solved)
%--------------------------------------------------------------------------
% I had a problem for passing a object tensor to the function. For some 
% reason when inside the function I can't evaluate M(2).
% This function coded within @tensor folder will be somthing like this:
% (not correct)
%
% % cast tensor-matrix into multidimensional array
% b = t.basis;
% sm = t.size;
% st = size(b, 1) * ones(1, t.order);
% nm = length(sm);
% nt = length(st);
% c = reshape(t.components, [sm st]);
% 
% % reorder dimensions and reshape to matrix
% c = permute(c, [1:nm/2 nm+(1:nt/2) nm/2+(1:nm/2) nm+nt/2+((1:nt/2))]);
% n = prod([sm(1:nm/2) st(1:nt/2)]);
% c = reshape(c, [n n]);
%--------------------------------------------------------------------------

% spatial dimension
dim  = size(b, 1);

% component order
order = o;

if order==0
S=t;
elseif order ==1 % TENSOR MATRIX WITH 1ST-ORDER TENSOR COMPONENTS
    % case 2D
    if dim==2
       
        % extract basis component
        e1 = b(1);
        e2 = b(2);
    
        % size of square tensor matrix nxn
        n   = size(t,1);
        % size of square tensor matrix nxn
        m   = size(t,2);
    
    
        % size of square scalar matrix NxM
        N   = dim * size(t,1);
        M   =  1  * size(t,2);
        
        % define scalar matrix S equivalent to t
        S = zeros(N,M);
        
        % SUBMATRICES
            % matrix with FIRST components of 1st-order tensor
            t1 = dot(t,e1);
    
            % matrix with SECOND components of 1st-order tensor
            t2 = dot(t,e2);
        
        % INDICES
            % row indices of component 1
            ii1 = 1:2:N-1;
    
            % row indices of component 1
            ii2 = 2:2:N;
    
            % column indices
            jj  = 1:1:M;
        
        % ASSIGN SCALAR MATRIX
            S(ii1,jj) = t1;
            S(ii2,jj) = t2;
      
    % case 3D
    elseif dim==3
       error("to be implemented");
    end

elseif order == 2   % TENSOR MATRIX WITH 2ST-ORDER COMPONENTS
    
    % case 2D
    if dim==2
       
        % extract basis component
        e1 = b(1);
        e2 = b(2);
    
        % size of square tensor matrix nxn
        n   = size(t,1);
        % size of square tensor matrix nxn
        m   = size(t,2);
    
    
        % size of square scalar matrix NxM
        N   = dim * size(t,1);
        M   = dim * size(t,2);
        
        % % define scalar matrix S equivalent to t
        % S = zeros(N,M);

        % define scalar matrix S equivalent to t
        S = sparse(N,M);
        
        % SUBMATRICES
            % matrix with components 11 of 2st-order tensors
            t11 = dot(e1,t,e1);
    
            % matrix with components 12 of 2st-order tensors
            t12 = dot(e1,t,e2);

            % matrix with components 21 of 2st-order tensors
            t21 = dot(e2,t,e1);

            % matrix with components 11 of 2st-order tensors
            t22 = dot(e2,t,e2);
        
        % INDICES
            % row indices of component 1
            ii1 = 1:2:N-1;
    
            % row indices of component 1
            ii2 = 2:2:N;
    
            % column indices
            jj1 = 1:2:M-1;

            % column indices
            jj2 = 2:2:M;
        
        % ASSIGN SCALAR MATRIX
            S(ii1,jj1) = t11;
            S(ii1,jj2) = t12;
            S(ii2,jj1) = t21;
            S(ii2,jj2) = t22;
    
    % case 3D
    elseif dim==3
         error("to be implemented");
    end


elseif order == 3   % TENSOR MATRIX WITH 3rd-ORDER COMPONENTS
       % case 2D
    if dim==2

    if size(t) ~= [1 1]; error('to be implemented'); end

    for i=1:2
        for j=1:2
            for k=1:2
                    S(i,j,k)=dot( b(j) , dot(b(i) , t , b(k)) );   
            end
        end
    end
    
    return;

    else
        error("to be implemented");
    end

elseif order == 4   % TENSOR MATRIX WITH 4th-ORDER COMPONENTS
    % case 2D
    if dim==2

    if size(t) ~= [1 1]; error('sorry, not supported for a matrix of 4th-order tensors'); end

    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    S(i,j,k,l)=dot( b(j) , dot(b(i) , t , b(l)), b(k));
                end
            end
        end
    end
    
    return;

    else
        error("to be implemented");
    end

else
    error("the tensor matrix component order is too high. tm2sm " + ...
          "function only transforms up to components of order 2");
end

% return a sparse matrix if it's bigger than 1x1
if abs(size(S))>3
S=sparse(S);
end

end