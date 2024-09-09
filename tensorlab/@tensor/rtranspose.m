function t2 = rtranspose(t1)

% RTRANSPOSE  Right-transponation of a tensor
%   B = RTRANSPOSE(A).
%   Return the right-transpose of a tensor. By definition, all the right
%   indices of the tensor are swapped. Therefore, only even order tensor
%   can be right-transposed.
%
%   For a matrix or multidimensional array A which contains tensors,
%   RTRANSPOSE(A) implies transposition of all tensor dimensions but not
%   the matrix/array dimensions.
%
%   See also TRANSPOSE, CTRANSPOSE, LTRANSPOSE.

if mod(t1.order,2) ~= 0
    error('Right transpose is only defined for an even order tensor')
end

% store components in multidimensional array
sm = t1.size;
st = size(t1.basis, 1)*ones(1, t1.order);
c1 = reshape(t1.components, [sm st]);

% permute array dimensions
nm = length(sm);
nt = length(st);
c2 = permute(c1, [1:nm nm+[1:1:nt/2 nt:-1:nt/2+1]]);

% store data in new tensor object
t2 = t1;
% t2.size = fliplr(sm);
t2.components = c2(:);

end