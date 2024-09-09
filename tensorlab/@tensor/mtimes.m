function t = mtimes(t1, t2)

% MTIMES   Dyadic or tensor product, use "*".
%   A*B, where A and B are vectors, determines the dyadic product, or
%   tensor product, AB and thus results in a second-order tensor.
%
%   A*B, where A and B may each be scalar, vectorial or tensorial,
%   determines their dyadic product. No contraction is performed and the
%   result is an nth-order tensor, where n is the sum of the number of
%   tensor dimensions of A and B.
%
%   A*B, where A and B are matrices of scalars, vectors or tensors,
%   returns the dyadic product of these matrices. The usual contraction
%   applies to the matrix dimensions of A and B and their inner dimensions
%   must correspond. If A or B is a 1x1 matrix, it is multiplied with all
%   components of the remaining matrix.
%
%   A*B*C*... returns the combined dyadic product of tensor-matrices A, B,
%   C, etc., i.e. ABC...
%
%   See also DYADIC, DOT, DDOT, CROSS.


% call function dyadic to compute dyadic product 
t = dyadic(t1, t2);

end