function t = cat(d, varargin)

t1 = varargin{1};
nt = size(t1.basis, 1)^t1.order;

for i = 1:nargin-1
    ti = varargin{i};
    c{i} = reshape(ti.components, [ti.size nt]);
end

c = cat(d, c{:});
s = size(c);
s = s(1:end-1);
if length(s) < 2
    s = [s 1];
end
c = c(:);

% store data in new tensor object
t = t1;
t.size = s;
t.components = c;
