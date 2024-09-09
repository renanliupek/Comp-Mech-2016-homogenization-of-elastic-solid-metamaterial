function t = sm2tm(S,b,o)
%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% TM2SM  Transform a an equivalent scalar matrix to tensor matrix
%   t = SM2TM(t,b,o)
%   
%   input  t - tensor matrix
%   input  b - basis
%   input  o - order of the tensor matrix component
%   output S - equivalent scalar matrix S 
%
%   ATTENTION: for
%--------------------------------------------------------------------------

% spatial dimension
dim  = size(b, 1);

% component order
order = o;

if order==0
t=S;

elseif order ==1 % SCALAR MATRIX WITH 1ST-ORDER TENSOR COMPONENTS
%     warning(['If the order==1, make sure that the first' ...
%         ' matrix index N of (N,M) will contract to a vector. For' ...
%         ' that, you' ...
%         ' might have to' ...
%         ' transpose the input of sm2tm().'])
    % case 2D
    if dim==2
       
        % extract basis component
        e1   = b(1);
        e2   = b(2);
    
        % size of eq. scalar matrix NxM
        N    = size(S,1);
        % size of eq. scalar matrix NxM
        M    = size(S,2);
    
    
        % size of square scalar matrix NxM
        n    = N/dim;
        m    = M ;
        
        % declare tensor matrix t 
        t    = zeros(o, b, n, m);
  
        
        % ASSIGN SCALAR MATRIX
        t    = S(1:2:N-1, :)*e1 + S(2:2:N, :)*e2;
      
    % case 3D
    elseif dim==3
       error("to be implemented");
    end

elseif order == 2   % SCALAR MATRIX WITH 2ST-ORDER COMPONENTS
    
    % case 2D
    if dim==2
       
        % extract basis component
        e1   = b(1);
        e2   = b(2);
    
        % size of eq. scalar matrix NxM
        N    = size(S,1);
        % size of eq. scalar matrix NxM
        M    = size(S,2);
               
    
        % size of square scalar matrix NxM
        n    = N/dim;
        m    = M/dim ;
        
        % declare tensor matrix t 
        t    = zeros(o, b, n, m);
  
        
        % ASSIGN SCALAR MATRIX
        t    = S(1:2:N-1, 1:2:M-1)*e1*e1 + S(1:2:N-1, 2:2:M)*e1*e2  + ...
               S(2:2:N  , 1:2:M-1)*e2*e1 + S(2:2:N  , 2:2:M)*e2*e2 ;

    
    % case 3D
    elseif dim==3
         error("to be implemented");
    end

else
    error("the tensor matrix component order is too high. tm2sm " + ...
          "function only transforms up to components of order 2");
end

end