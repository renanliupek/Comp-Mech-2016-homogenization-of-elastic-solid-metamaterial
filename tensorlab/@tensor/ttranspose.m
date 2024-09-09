function t2 = ttranspose(t1)

% TTRANSPOSE  Transpose of a tensor
%   B = TTRANSPOSE(A)
%     Gives the tensor transpose a A. For a single tensor the result is equal
%     to A' or A.'. For a matrix however, only the tensors are transposed.
%     The matrix is leaved untouched.
%
%     See also TRANSPOSE, CTRANSPOSE, RTRANSPOSE.

% store components in multidimensional array
sm = t1.size;
st = size(t1.basis, 1)*ones(1, t1.order);
c1 = reshape(t1.components, [sm st]);

% permute array dimensions
nm = length(sm);
nt = length(st);
c2 = permute(c1, [1:nm nm+(nt:-1:1)]);

% store data in new tensor object
t2 = t1;
% t2.size = fliplr(sm);
t2.components = c2(:);

end