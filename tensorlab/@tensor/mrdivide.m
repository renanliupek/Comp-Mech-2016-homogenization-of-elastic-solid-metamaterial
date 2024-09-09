function t = mrdivide(varargin)

% MRDIVIDE Vector/tensor dot division.
%   MRDIVIDE(A,B), where A and B are vectors, gives an error. It is
%   non-existent.
%
%   DYADIC(A,B), where A a vector/tensor and B is a scaler gives the same
%   result as A * (1/B)
%
%   DYADIC(A,B), where A a scalar and B is a vector/tensor gives an error
%   since it is not defined.
%
%
% call product function to compute dyadic division

for i = length(varargin)-1
    ti  = varargin{i};
    ti1 = varargin{i+1};
    if isa(ti,'tensor')
        if isa(ti1,'tensor')
            error('The dyadic division between two tensors is not defined')
        end
    else
        if isa(ti1,'tensor')
            error('A scalar divided by a tensor is not defined')
        end
    end
end

varargin{2} = 1 / varargin{2};

t = product(0, varargin{:});

end