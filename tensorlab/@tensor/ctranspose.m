function t2 = ctranspose(t1)

% CTRANSPOSE Tensor transposition.
%   A' is the tensor transpose of the tensor A.
%
%   For a matrix or multidimensional array A which contains tensors, A'
%   implies transposition of all tensor dimensions as well as the
%   matrix/array dimensions.
%
%   See also TRANSPOSE, LTRANSPOSE, RTRANSPOSE.

% store components in multidimensional array
sm = t1.size;
st = size(t1.basis, 1)*ones(1, t1.order);
c1 = reshape(t1.components, [sm st]);

% permute array dimensions
nm = length(sm);
nt = length(st);
c2 = permute(c1, [nm:-1:1 nm+(nt:-1:1)]);

% store data in new tensor object
t2 = t1;
t2.size = fliplr(sm);
t2.components = c2(:);

end