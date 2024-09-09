function varargout = scatter(X,varargin)

% SCATTER Scatter/bubble plot
%    SCATTER(X) plots markers at the positions indicated by X.
%
%    SCATTER(X,S,U) plots colored markers at the positions X. The size of
%    the markers (in points^2) is determined by the column matrix S which
%    must comprise only positive numbers. The color, on the other hand,
%    is determined by the column U (which may comprise both positive
%    and negative values). The values in U are linearly mapped to the
%    colors in the current colormap. A legend for the different colors
%    can therefore be added by subsequently executing:
%      scatter(X,S,U) 
%      colorbar
%    notice how X, S, and U are all of the same matrix dimensions.
%
%    SCATTER(X,[],U) plots colored markers at the positions X which are
%    all of the same size.
%    
%    SCATTER(...,'filled') fills the markers.
%
%    All properties of the regular Matlab scatter are applicable for its
%    Tensorlab counterpart. 
%
%    SEE ALSO: SCATTER, TENSOR/QUIVER
% 

% Check if X has a tensorial character and is of order 1
if ~isa(X, 'tensor')
    error('X must be a Tensor')
elseif X.order ~= 1
    error('X must be or order one (i.e. a vector)')
end

x  = X.components;
xd = size(X.basis,1);
xn = size(x,1);
x  = reshape(x, xn/xd , xd );

if xd == 2
    
    h = scatter(x(:,1),x(:,2),varargin{:});
    
elseif xd == 3
    
    h = scatter3(x(:,1),x(:,2),x(:,3),varargin{:});
    
end

% set output
if     nargout==1, varargout{1}=h;
elseif nargout~=0, error('Too many output arguments requested'); 
end

end
