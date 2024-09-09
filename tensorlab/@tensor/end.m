function e = end(t, i, n)

% linear indexing
if n == 1
    e = prod(t.size);
    
% multiple indices
else
    e = t.size(i);
    
end
