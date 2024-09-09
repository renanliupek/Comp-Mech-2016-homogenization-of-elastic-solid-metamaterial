function new_indices = get_ind(conversion_rule, old_indices)
% Renan Liupekevicius TU/e
% GET_IND find the new indices of 'old_indices' given 'conversion_rule'
% 
%   new = GET_IND(conversion_rule, old_indices)
%   
%   conversion_rule - old indices are the elements and new indices are the
%                     corresponding vector index;
%
%   old_indices     - set o nodes to be converted to new index definition

if isempty(conversion_rule) new_indices=[]; 
    return;
    warning('returned empty index list');
end

new_indices = zeros(1,length(old_indices));

for i=1:length(old_indices)
    new_indices(i) =  find( conversion_rule == old_indices(i) );
end

end