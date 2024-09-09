function varargout = eig(t)

% EIG    Eigenvalues and eigenvectors of a tensor or matrix of tensors.
%    E = EIG(A), where A is a tensor object, returns a column matrix
%    containing the eigenvalues of A. The tensor A must be of even order.
%    If A is a matrix of tensors, the matrix must be square and the tensors
%    it contains must be of even order.
%
%    [V, D] = EIG(A) computes the eigenvalues and eigenvectors of A. The
%    diagonal matrix D has the eigenvalues on its main diagonal. The
%    columns of matrix V contain the corresponding eigenvectors so that
%    A*V = V*D. Note that, depending on the tensorial order of A, matrix V
%    may be scalar, vectorial or tensorial; its tensorial order equals that
%    of A divided by two. 

% check input argument
if ~isa(t, 'tensor')
    error('Argument is not a tensor.')
end
if mod(t.order, 2) > 0
    error('Tensor should be of even order.')
end
if length(t.size) > 2
    error('Matrix should be two-dimensional.')
end
if t.size(1) ~= t.size(2)
    error('Matrix must be square.')
end

% check number of output arguments
if nargout > 2
    error('Too many output arguments.')
end

% cast tensor-matrix into multidimensional array
b = t.basis;
sm = t.size;
st = size(b, 1) * ones(1, t.order);
nm = length(sm);
nt = length(st);
c = reshape(t.components, [sm st]);

% reorder dimensions and reshape to matrix
c = permute(c, [1:nm/2 nm+(1:nt/2) nm/2+(1:nm/2) nm+nt/2+((1:nt/2))]);
n = prod([sm(1:nm/2) st(1:nt/2)]);
c = reshape(c, [n n]);

% compute eigenvalues and, if necessary, eigenvectors
if nargout <= 1
    
    % compute eigenvalues and store them in output argument
    lambda = eig(c);
    varargout(1) = {lambda};
    
else
    
    % compute eigenvalues and components of eigenvectors
    [c, lambda] = eig(c);
    
    % reorder components to form multidimensional array
    c = reshape(c, [sm(1:nm/2) st(1:nt/2) n]);
    c = permute(c, [1:nm/2 nm/2+nt/2+1 nm/2+(1:nt/2)]);
    
    % cast eigenvectors into tensor object
    v = tensor;
    v.basis = b;
    v.size = [sm(1:nm/2) n];
    v.order = nt/2;
    v.components = c(:);

    % store both results in output arguments
    varargout(1) = {v};
    varargout(2) = {lambda};

end

end
