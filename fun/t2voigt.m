function Cv = t2voigt(C,b,o)

%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e
% T2VOIGT  Transform a tensor (array) to an scalar matrix using Voigt 
% notation
%
%   M = T2VOIGT(t,b,o)
%   
%   input  C  - tensor
%   input  b  - basis
%   input  o  - order of the tensor
%   output Cv - scalar matrix  
%
%--------------------------------------------------------------------------

% spatial dimension
dim  = size(b, 1);

if o==4


% case 2D
if dim==2

    % define indices 
    disp('index pairs on voigt notation')
    index_pair={[1 1], [2 2], [2 1], [1 2]};
    
    % initialize output voigt matrix
    Cv = zeros(length(index_pair),length(index_pair));

    for IJ = 1:length(index_pair)
        for KL = 1:length(index_pair)
            Cv(IJ,KL)= C(index_pair{IJ}(1),index_pair{IJ}(2),...
                         index_pair{KL}(1),index_pair{KL}(2));
        end
    end




% case 3D
elseif dim==3
    

    % define indice pair
    disp('index pairs on voigt notation');
    index_pair={[1 1], [2 2], [3 3], [2 1], [1 2], [3 1], [1 3], [3 2], [2 3]}
    
    % initialize output voigt matrix
    Cv = zeros(length(index_pair),length(index_pair));


    for IJ = 1:length(index_pair)
        for KL = 1:length(index_pair)
            Cv(IJ,KL)= C(index_pair{IJ}(1),index_pair{IJ}(2),...
                                          index_pair{KL}(1),index_pair{KL}(2));
        end
    end
        

else 
    error('t2voigt supports only 2D or 3D tensors');
end



end

end
