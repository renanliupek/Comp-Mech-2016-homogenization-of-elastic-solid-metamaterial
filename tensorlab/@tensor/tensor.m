function t = tensor(s)

if nargin == 0
    
    t.basis = '';
    t.order = 0;
    t.size = [0 0];
    t.components = [];
    t = class(t, 'tensor');
    
elseif isa(s, 'tensor')

    t = s;
    
elseif isstruct(s)
    
    t = class(s, 'tensor');
    
else
    
    error('Illegal input argument')
    
end
