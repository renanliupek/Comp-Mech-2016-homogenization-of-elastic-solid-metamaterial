function t2 = inv(t1)

% INV  Inverse of a tensor.
% 
%   INV(C) caculates the inverse of a tensor object. It is only possible to
%   calculate the inverse of first order tensors (i.e. vectors) or second
%   order tensors.
% 

% Check the input
if t1.order > 2
    error('The inverse is only defined for first and second order tensors');
end

% Compute the inverse for second order tensor. First make a matrix out if
% it, calculate the inverse, and then make a tensor again.
if t1.order == 2
    
    st = size(t1.basis, 1)*ones(1, t1.order);
    c = reshape(t1.components, st);
    c = inv(c);
    c = reshape(c,size(c,1)*size(c,2),1);
    t2 = t1;
    t2.components = c;
    
% Compute the inverse for first order tensors (vectors). The inverse is now
% simply the inverse of each component.
elseif t1.order == 1 
    
    t2 = t1;
    
    for ii = 1:size(t2.components,1)
        t2.components(ii) = 1/t2.components(ii);
    end
    
end

% c2 = t1.components;
% c2 = [c2(1) c2(3);
%       c2(2) c2(4)];
% c2 = inv(c2);
% c2 = [c2(1,1); c2(2,1); c2(1,2); c2(2,2)];
% t2 = t1;
% t2.size = t1.size;
% t2.components = c2(:);