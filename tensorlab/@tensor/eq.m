function [ bool ] = eq ( A , B )


if ~isa(A,'tensor') || ~isa(B,'tensor')
%     warning('Both inputs must be a tensor object')
    bool = 0;
elseif A.order ~= B.order
%     warning('Inputs have a different order')
    bool = 0; 
elseif A.basis ~= B.basis
%     warning('Inputs have a different basis')
    bool = 0;
elseif A.size ~= B.size
%     warning('Inpunt have a different size')
    bool = 0;
else
    bool = 1;
    for i = 1:size(A.components,1)
        if A.components(i) ~= B.components(i)
            bool = 0;
        end
    end
end

end

