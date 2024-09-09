function c = det(t1)

% DET  Determinant of a tensor.
% 
%   DET(C) caculates the inverse of a tensor object. It is only possible to
%   calculate the determinant of a second order tensor.
% 

% Check the input
if t1.order ~= 2
    error('The determinant is only defined for second order tensors');
end

% Compute the determinant for second order tensor. First make a matrix out
% if it. Then the determinant of this matrix is caluculated
    
st = size(t1.basis, 1)*ones(1, t1.order);
c = reshape(t1.components, st);
c = det(c);


end

% c1 = t1.components;
% c1 = [c1(1) c1(3);
%       c1(2) c1(4)];
% c2 = det(c1);