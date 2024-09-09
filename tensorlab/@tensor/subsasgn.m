function t3 = subsasgn(t1, s, t2)

% check type of subcript
if s.type ~= '()'
    error('Illegal type of subscript.')
end

% check consistency of input matrices
if ~all(t1.basis(:) == t2.basis(:))
    error('Tensors in subscripted assignment should have same basis.')
elseif t1.order ~= t2.order
    error('Tensors in subscripted assignment should be of the same order.')
end

% check number of indices
ss = s.subs;
if ~( ((length(ss) == 1) && (length(t2.size) == 2)) || ...
        ((length(ss) == length(t1.size)) && (length(ss) == length(t1.size))) )
    error('Number of indices does not match matrix dimensions.')
end

% replace components
if length(ss) == 1;
    c1 = reshape(t1.components, [prod(t1.size) size(t1.basis, 1)^t1.order]);
    c2 = reshape(t2.components, [prod(t2.size) size(t2.basis, 1)^t2.order]);
else
    c1 = reshape(t1.components, [t1.size size(t1.basis, 1)^t1.order]);
    c2 = reshape(t2.components, [t2.size size(t2.basis, 1)^t2.order]);
end
ss = [ss ':'];
c1(ss{:}) = c2;
n = size(c1);
n = n(1:end-1);
if length(n) < 2
    n = [n 1];
end
c1 = c1(:);

% store data in new tensor object
t3 = t1;
t3.size = n;
t3.components = c1;
