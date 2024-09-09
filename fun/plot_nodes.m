function plot_nodes(x,nodes,ee)
%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% PLOT_NODES  plot nodes if tensor vector x using basis ee
% 
%   PLOT_NODES(t,b,o)
%   
%   input x     - tensor vector
%   input nodes - array of node numbers to be plotted according to x index
%   input ee    - basis
%   output      - a graphic
%
%--------------------------------------------------------------------------
plot(dot(x(nodes),ee(1)),dot(x(nodes),ee(2)),'*')

end