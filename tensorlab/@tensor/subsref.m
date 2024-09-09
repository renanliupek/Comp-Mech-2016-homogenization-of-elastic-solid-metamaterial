function t2 = subsref(t1, s)

% check type of subcript
if s.type ~= '()'
    error('Illegal type of subscript.')
end

% check number of indices
ss = s.subs;
if (length(ss) ~= 1) && (length(ss) ~= length(t1.size))
    error('Number of indices does not match matrix dimensions.')
end

% extract components
if length(ss) == 1;
    c = reshape(t1.components, [prod(t1.size) size(t1.basis, 1)^t1.order]);
else
    c = reshape(t1.components, [t1.size size(t1.basis, 1)^t1.order]);
end
ss = [ss ':'];
c = c(ss{:});
n = size(c);
n = n(1:end-1);
if length(n) < 2
    n = [n 1];
end
c = c(:);

% store data in new tensor object
t2 = t1;
t2.size = n;
t2.components = c;

end
