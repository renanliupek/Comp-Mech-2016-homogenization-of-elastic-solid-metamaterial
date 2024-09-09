function plot(X,varargin)

% PLOT(X,U)  Plot of tensor objects
% 
%   PLOT(X,U,Properties)
%     Plots (a column) U versus (a column) X. X has to be a (column of) 
%     vector(s). Depending on the order and basis of U, the resulting plot 
%     is different. PLOT(X,U) calls the correct plotting function.
%     
%     If U is a column of scalars, the resulting plot is a scatter plot, 
%     where the size and color of the dots indicate the corresponding
%     value. Plot uses the scatter function with some options already
%     specified. The function 
%       scatter(X,S,C,'filled')
%     is called. Here, S = C = U. Other options can still be specified.
%     Consult <a href = "matlab:help tensor/scatter">help tensor/scatter</a> for more information.
%
%     If U is a column of vectors, the resulting plot is a quiver plot,
%     where the vectors of U are plotted at starting positons X. If no
%     color is provided, the color is generated automatically, such that
%     the it represents the length of the vector.  Plot uses the quiver
%     function such that the scaling is off. I.e., the function
%       quiver(X,U,0)
%     is called. Other options can still be specified. 
%     Consult <a href = "matlab:help tensor/quiver">help tensor/quiver</a> for more information.
%
%     If U is a column of 2nd order tensors, the resulting plot is a quiver
%     plot, where the vectors correspond to the eigenvectors of the
%     components of U, with length equal to their eigenvalues. If no
%     color is provided, the color is generated automatically, such that
%     the it represents the respective eigenvalue.  Plot uses the quiver
%     function such that the scaling is off. I.e., the function
%       quiver(X,U,0)
%     is called. Other options can still be specified. 
%     Consult <a href = "matlab:help tensor/quiver">help tensor/quiver</a> for more information.
%                                                            

warning('tensor:new','Please use scatter of quiver directly,\nthis function will be removed in newer functions of tensorlab');

if nargin > 1
    
    U = varargin{1};
    varargin = { varargin{ 2:size(varargin,2) } };
    
    if isa(U,'tensor')
        
        quiver(X,U,0,varargin{:});
        
    else
        
        Smin = 1;
        Smax = 2;
        S = (max(U)-U)/(max(U)-min(U)) * Smin + ...
            (U-min(U))/(max(U)-min(U)) * Smax;
        S = 40*S;
        
        scatter(X,S,U,varargin{:},'filled');
        
        colorbar
        
    end
    
else
    
    scatter(X,varargin{:});
    
end
