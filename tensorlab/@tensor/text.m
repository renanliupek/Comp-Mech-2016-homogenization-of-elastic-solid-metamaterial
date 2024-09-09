function text ( x , name , varargin )

% TEXT  Text annotation.
%   TEXT(COORD,NAME), where COORD is a vector (of the tensor class) that
%   gives the position of the text in the figure. NAME is a string
%   containing the text label you wish to assign.
% 
%   TEXT(COORD,NAME,PROPERTIES), where PROPERTIES are the same properties
%   you can specify in the regular MATLAB text function (type help text, to
%   get more information).
% 
%   Consult <a href = "matlab:help text">help text</a> for information on the regular Matlab text.

if ~isa(x, 'tensor')

    error('Coord should have a tensorial character.')
    
elseif length(x.basis) == 1
    
    error('Coord should have a 2D of 3D basis.')
     
end

basis = x.basis;
dim   = size(basis,1);
x     = x.components;
nx    = size(x,1)/dim;
x     = reshape(x,nx,dim);

text('Position',x,'String',name,varargin{:});

end
