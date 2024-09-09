function t = double2tensor(d, e)

s.basis = e;
s.order = 0;
s.size = size(d);
s.components = d(:);

t = tensor(s);

