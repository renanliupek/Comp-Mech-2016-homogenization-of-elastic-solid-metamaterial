function h=patch ( varargin )
% 
% PATCH  Create Patch.
%   PATCH( ... ) is the tensor equivalent of the default Matlab patch. The
%   function constructs planes. The difference with the default Matlab
%   patch is that for the tensor equivalent the option 'Vertices' has to be
%   specified explicitly. Moreover, the 'Vertices' entries have to be
%   vectors (i.e. position vectors).
% 
% Consult <a href = "matlab:help patch">help patch</a> for more information on the regular Matlab patch.
% 


% Find the vertices, which will be containin the position vectors
i = 1;
while ~strcmp(varargin{i},'Vertices') && i < nargin
    i = i + 1;
end

% Reshape the column of vectors to have \vec{x}_i = [x_i y_i z_i]
if strcmp(varargin{i},'Vertices') && isa(varargin{i+1},'tensor')
        
    x   = varargin{i+1};
    dim = size(x.basis,1);
    x   = x.components;
    nx  = size(x,1)/dim;
    x   = reshape(x,nx,dim);
    
    % Group all input options, except 'Vertices',x
    if i+2 > nargin
        options = { varargin{1:i-1} };
        
    else
        options = { varargin{1,1:i-1},varargin{1,i+2:nargin} };
        
    end
    
    % Make the actual plot
    h=patch('Vertices',x,options{:});
    
else
    
    % If no tensor object was encountered
    h=patch(varargin);
    
end

end
    
    



