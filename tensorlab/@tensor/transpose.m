function t2 = transpose(t1)

% TRANSPOSE  Matrix transposition of tensor-matrix.
%   A.' or TRANSPOSE(A) is the matrix transpose of the tensor-matrix A. No
%   transposition is applied to the tensor-elements of A. This is opposed
%   to A' or CTRANSPOSE(A), which does imply transposition of the tensors
%   as well.
%
%   See also CTRANSPOSE.

% store components in multidimensional array
sm = t1.size;
st = size(t1.basis, 1)*ones(1, t1.order);
c1 = reshape(t1.components, [sm st]);

% permute array dimensions
nm = length(sm);
nt = length(st);
c2 = permute(c1, [nm:-1:1 nm+(1:nt)]);

% store data in new tensor object
t2 = t1;
t2.size = fliplr(sm);
t2.components = c2(:);

end