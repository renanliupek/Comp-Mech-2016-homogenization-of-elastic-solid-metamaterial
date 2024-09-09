function [ H ] = surf ( varargin )

if isfloat(varargin{1})
    
    h = varargin{1};
    X = varargin{2};
    
else
    
    X = varargin{1};
    
end

if ~isa(X,'tensor')
    error('The coordinates must be tensor objects\n');    
end


x  = X.components;
xd = size(X.basis,1);
xn = size(x,1);
x  = reshape(x, xn/xd , xd );
    
if exist('h') ~= 1
    
    varargin{2:end}
    H = surf(x(:,1),x(:,2),varargin{2:end});
    
else

    H = surf(h,x(:,1),x(:,2),x(:,3),varargin{3:end});
    
end


end