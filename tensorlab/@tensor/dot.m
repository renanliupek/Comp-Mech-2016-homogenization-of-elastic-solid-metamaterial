function t = dot(varargin)

% DOT    Vector/tensor dot product.
%   DOT(A,B), where A and B are vectors, computes the dot
%   product, or inner product, A.B and thus results in a scalar value.
%
%   DOT(A,B), where A and B are tensors, computes the dot product of these
%   tensors, i.e. a single contraction is performed between the final
%   tensor dimension of A and the first of B.
%
%   DOT(A,B), where A and B are matrices of vectors or tensors,
%   returns the dot product of these matrices. The usual contraction
%   applies to the matrix dimensions of A and B and their inner dimensions
%   must correspond. If A or B is a 1x1 matrix, it is multiplied with all
%   components of the remaining matrix.
%
%   DOT(A,B,C,...) returns the combined dot product of tensor-matrices A,
%   B, C, etc., i.e. A.B.C. ...
%
%   See also DYADIC, CROSS, DDOT.

% call product function to compute dot product
t = product(1, varargin{:});

end