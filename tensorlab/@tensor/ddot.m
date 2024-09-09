function t = ddot(varargin)

% DDOT   Vector/tensor double dot product.
%   DDOT(A,B), where A and B are tensors of at least order two, computes
%   the double dot product, or double inner product, A:B. A double
%   contraction is performed between the penultimate and ultimate tensor
%   dimesions of A and respectively the first and second tensor dimensions
%   of B. Depending on the tensor characters of A and B, a scalar, vector
%   or tensor is returned.
%
%   DDOT(A,B), where A and B are matrices of vectors or tensors,
%   returns the double dot product of these matrices. The usual contraction
%   applies to the matrix dimensions of A and B and their inner dimensions
%   must correspond. If A or B is a 1x1 matrix, it is multiplied with all
%   components of the remaining matrix.
%
%   DDOT(A,B,C,...) returns the combined double dot product of
%   tensor-matrices A, B, C, etc., i.e. A:B:C. ...
%
%   See also DYADIC, DOT, CROSS, PRODUCT.

% call product function to compute double dot product
t = product(2, varargin{:});

end