function v=treal(v,ee)
%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% TREAL  real part of a 1st-order tensor (2D)
%   v = TREAL(t,b)
%   
%   input  t - tensor 
%   input  b - basis
%   output v - real part of a tensor
%
%--------------------------------------------------------------------------

v = real(dot(v,ee(1)))*ee(1)+ real(dot(v,ee(2)))*ee(2);

end