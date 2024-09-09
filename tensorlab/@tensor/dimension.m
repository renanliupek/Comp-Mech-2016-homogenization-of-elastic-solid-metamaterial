function [ dim ] = dimension ( varargin )

% DIMENSION  Dimension of a tensor.
% 
%     [dim] = DIMENSION ( X )
%      This function returns the dimension of a tensor object. I.e. it return
%      whether a tensorobject (e.g. a vector) has a two-dimension or
%      three-dimensional basis. 
% 
%     [dim] = DIMENSION ( X1 , X2 , ... )
%      Returns the dimension that all input tensorobjects share. An error is
%      given when the dimension is not the same for all objects. 


for i = 1:size(varargin,2)
    if i == 1
        dim = size(varargin{i}.basis,1);
    elseif size(varargin{i}.basis,1) ~= dim;
        error('Not all input values have the same dimension')
    end
end