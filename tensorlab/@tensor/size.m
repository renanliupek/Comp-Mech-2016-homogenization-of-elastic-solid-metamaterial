function s = size(t, d)

if nargin < 2

    s = t.size;

else

    if 1 <= d <= length(t.size)
        s = t.size(d);
    else
        error('Index exceeds matrix dimensions.')
    end
    
end