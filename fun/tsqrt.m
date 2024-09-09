function v=tsqrt(v,ee)
%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% TSQRT  square root of a 1st-order tensor element-wise (2D)
%   v = TSQRT(t,b)
%   
%   input  t - tensor 
%   input  b - basis
%   output v - sqrt of a tensor
%
%--------------------------------------------------------------------------

v = sqrt(dot(v,ee(1)))*ee(1)+ sqrt(dot(v,ee(2)))*ee(2);

end