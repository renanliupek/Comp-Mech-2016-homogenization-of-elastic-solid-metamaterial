function t = minus(t1, t2)

if ~isa(t1, 'tensor') || ~isa(t2, 'tensor')

    error('Tensorial characters must agree.')

elseif size(t1.basis, 1) ~= size(t2.basis, 1) || any(any(t1.basis ~= t2.basis))

    error('Tensor bases must agree.')

elseif t1.order ~= t2.order

    error('Tensiorial orders must agree')

elseif any(t1.size ~= t2.size)

    error('Matrix dimensions must agree')

else

    t = t1;
    t.components = t1.components - t2.components;
    
end
