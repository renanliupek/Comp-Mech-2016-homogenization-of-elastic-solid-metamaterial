function [ true ] = checkbasis ( e )

% CHECKBASIS  Check if basis is orthogonal.
%      DIM = CHECKBASIS ( e )
%      checks the dimension of the basis, and wheter of not is is orthogonal.
%      The output is the dimension of the basis.
% 
    
if ~isa( e ,'tensor')
    error('The basis must be a tensor')
    
elseif size(e,1) > 3
    error('The maximum dimension of the basis is 3')
    
elseif size(e,1) == 1
    if size(e.basis,1) == 1
        true = 1;
    else
        error('The number of basis components does not match the basis')
    end

elseif size(e,1) == 2
    if size(e.basis,1) == 2

        S.type = '()';
        S.subs = {[1]};
        e1 = subsref(e,S);
        S.subs = {[2]};
        e2 = subsref(e,S);

        if dot(e1,e2) ~= 0
            error('The basis is not othogononal')
        end
    
        true = 2;
        
    else
        error('The number of basis components does not match the basis')
    end

elseif size(e,1) == 3
    if size(e.basis,1) == 3
        
        S.type = '()';
        S.subs = {[1]};
        e1 = subsref(e,S);
        S.subs = {[2]};
        e2 = subsref(e,S);
        S.subs = {[3]};
        e3 = subsref(e,S);

        if dot(e1,e2) ~= 0 || dot(e1,e3) ~= 0 || dot(e2,e3) ~= 0
            error('The basis is not othogononal') 
            true = 0;
        end

        true = 3;
        
    else
        error('The number of basis components does not match the basis')
    end

end
       
    
end