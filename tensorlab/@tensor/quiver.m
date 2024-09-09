function varargout = quiver(X,U,varargin)
% QUIVER Quiver plot.
%
%    QUIVER(X,U) quiver plot. If 
% 
%     - U is (a column of) vectors: the components of U are plotted as
%       vectors, starting at the position(s) provided by X. QUIVER
%       automatically generate colors such that they correspond to length
%       of the vectors U.
% 
%     - U is (a column of) second order tensors: the eigenvectors of U
%       are plotted. The length of these vectors is determined by the
%       eigenvalues. QUIVER automatically generate colors such that they
%       correspond to eigenvalues of U.
% 
%    QUIVER(X,U,0) omits the automatic scaling of QUIVER.
% 
%    All regular Matlab quiver properties are applicable for its Tensorlab
%    counterpart.
%
% SEE ALSO: QUIVER, TENSOR/SCATTER

% Checking the input

% Check if the basis is 2D or 3D
if (size(U.basis,1) ~= 2) && (size(U.basis,1) ~= 3)
    error('Plot of 1D vector undefined')
end

% Check if X has a tensorial character and is of order 1
if ~isa(X, 'tensor')
    error('X must be of a tensorial character')
elseif X.order ~= 1
    error('X must have a vectorial character')
end

% Check if U has the same dimensions as X
if size(U) ~= size(X)
    error('Matrix dimensions of U and X do not agree')
end

% Check if U is of order 2 or less and if the bases of X and U are the same
if isa(U, 'tensor')
    if (U.order > 2)
        error('U must have order 2 or less')
    elseif size(U.basis, 1) ~= size(X.basis, 1) || any(any(U.basis ~= X.basis))
        error('Tensor bases must agree')
    end
end

if size(X,2) > 1
    if size(X,1) > 1
        error('X must be an array')
    end
    X = X.';
end

% General parameters
dim = size(X.basis,1);
sup = 0;

% Find whether or not a color is specified. If one is specified, remove in
% from the varargin.
nin = size(varargin,2);
ii = 1;
while ii < nin && ...
      ~( strcmp(varargin{ii},'Color') || strcmp(varargin{ii},'color') )
  ii = ii + 1;
end

if ii > 0 && ii < nin
    if strcmp(varargin{ii},'Color') || strcmp(varargin{ii},'color')
        sup = 1;
        ii  = ii + 1;
        c   = varargin{ii};
        varargin = { varargin{ [1:ii-2 ii+1:nin] } };
        if length(c) ~= size(X,1)
            cm = {};
            if U.order == 1
                list = 1:size(U,1);
            elseif U.order == 2
                list = 1:size(U,1)*dim;
            end

            for jj = list
                cm = {cm{:} c}; % Make a cell with the constant color for each vector
            end
        else
            cm = c;
            
            if ~isa(cm,'cell')
                error('Unknown color input')
            end
        end
    end
elseif ii == nin && ( strcmp(varargin{ii},'Color') || strcmp(varargin{ii},'color') )
    error('No color specified');
else
    cm = cmap(U,'jet');
end

x  = X.components;
xd = size(X.basis,1);
xn = size(x,1);
x  = reshape(x, xn/xd , xd );

if U.order == 1
    u  = U.components;
    ud = size(U.basis,1);
    un = size(u,1);
    u  = reshape(u, un/ud , ud );
    
elseif U.order == 2
    
    u = zeros(dim*size(U,1),dim);
    
    S.type = '()';
    jj = 1;
    for ii = 1:size(U,1)
        S.subs = {[ii]};
        Ui = subsref(U,S);
        [Ui lambda] = eig(Ui);
        lambda = diag(lambda);
        
        ui  = Ui.components;
        uid = size(Ui.basis,1);
        uin = size(ui,1);
        ui  = reshape(ui, uin/uid , uid );

        for i = 1:dim
            u(jj,:) = ui(i,:)*lambda(i); jj = jj + 1;
        end
    end     
    
    xt = x;
    if size(U.basis,1) == 2
        x = zeros(size(xt,1)*2,size(xt,2));
        jj = 1;
        for ii = 1:size(xt,1)
            x(jj,:) = xt(ii,:); jj = jj + 1;
            x(jj,:) = xt(ii,:); jj = jj + 1;
        end
    elseif size(U.basis,1) == 3
        x = zeros(size(xt,1)*3,size(xt,2));
        jj = 1;
        for ii = 1:size(xt,1)
            x(jj,:) = xt(ii,:); jj = jj + 1;
            x(jj,:) = xt(ii,:); jj = jj + 1;
            x(jj,:) = xt(ii,:); jj = jj + 1;
        end
    end
    
end


hold on

Handle = zeros(size(x,1),1);

if size(U.basis,1) == 2
    
    for ii = 1:size(x,1)
              
        Handle(ii) = ...
         quiver(x(ii,1),x(ii,2),u(ii,1),u(ii,2),varargin{:},'Color',cm{ii});
        
    end
    
elseif size(U.basis,1) == 3
    
    for ii = 1:size(x,1)
              
        Handle(ii) = ...
         quiver3(x(ii,1),x(ii,2),x(ii,3),u(ii,1),u(ii,2),u(ii,3),varargin{:},'Color',cm{ii});
        
    end
    
end

hold off
if ~sup
    colorbar
end

% set output
if     nargout==1, varargout{1}=Handle;
elseif nargout~=0, error('Too many output arguments requested'); 
end

end
