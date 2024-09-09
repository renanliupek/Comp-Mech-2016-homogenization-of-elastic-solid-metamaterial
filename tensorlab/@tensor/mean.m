function [ xbar ] = mean ( coord )
% 
% MEAN  The mean of a tensor (or vector).
% 
%     xbar = MEAN ( x ) calculates the average of a (or multiple) column(s) of
%     vectors, a (or multiple) column of tensors
% 
%     If x is a column the output will simply be the average of this column.
%     The components of the avarge vector or tensor (depending on the input)
%     are the average of the components of all vectors (or tensors) in the
%     column x
% 
%     If x are multiple columns of vectors (or tensors). The output is now a row
%     vector where each enty is the average (same as previous) of that column.
%     E.g. xbar(1,3) = mean( x(:,3) ) etc.
% 

% % Check the input
S.type = '()';
S.subs = {[1]};
coord1 = subsref(coord,S);

if isa(coord, 'tensor')
    if size(coord,1) > 1
        
        for i = 1:size(coord,1)-1
            
            S.type = '()';
            S.subs = {[i]};
            coordi = subsref(coord,S);
            S.type = '()';
            S.subs = {[i+1]};
            coordiplus = subsref(coord,S);

            if size(coordiplus.basis, 1) ~= size(coordi.basis, 1) ||...
               any(any(coordiplus.basis ~= coord1.basis))
                if coord.order == 1
                    error('Vector bases must agree.')
                else
                    error('Tensor bases must agree.')
                end
            end

            if coordiplus.order ~= coord1.order
                error('All input must be of the same order')
            end

        end
    end
        
        
elseif ~isnumeric(coord)
    error('Mean of tensor and %s is undefined.', class(coord))
end
if ndims(coord) > 2
    error('Input arguments must be two-dimensional matrices.')
end

if isempty(coord1.basis)
    error('No tensor object found in product.')
end

% Calculate the mean

basis  = coord.basis;
dim    = size(basis,1);
order  = coord.order;
x      = coord.components;
orsize = coord.size;

x      = reshape(x,dim*orsize(1)*order,orsize(2));
xmean  = zeros(dim*order,orsize(2));

for ww = 1:orsize(2)
    j = 1;
    for i = 1:length(xmean)
        xmean(i,ww) = mean(x([j:j+orsize(1)-1],ww));
        j = j+orsize(1);
    end
end

coord.size(1) = 1;
coord.components = zeros(dim*order*orsize(2),1);

for i = 1:orsize(2)
    coord.components(1+(i-1)*dim*order:i*dim*order,1) = xmean(:,i);
end

xbar = coord;

end