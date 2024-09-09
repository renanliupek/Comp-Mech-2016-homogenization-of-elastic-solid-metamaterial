function [ s ] = norm ( x , varargin )

% NORM  Vector norm.
%   NORM( X )
%   Return the norm of vector X. In this case the Eucledian 2-norm, or
%   length of the vector is calculated. X is restricted to a vector (i.e.
%   the maximum order is 1).
% 
%   NORM( X , t )
%   t specifies which norm to use. The choices are
%     t = 1    Eucledian 1-norm (sum of absolute values of components)
%     t = 2    Eucledian 2-norm (length of vector)
%     t = inf  maximum absolute value of the components
% 

if size(varargin,2) == 0
    n = 2;
elseif size(varargin,1) ~= 1 || size(varargin,2) ~= 1
    error('Unknown input')
else
    n = varargin{1};
    if n ~= 1 && n ~= 2 && n ~= inf
        error('Only the 1, 2, or inf norm are defined')
    end
end

if ~isa(x,'tensor')
    error('first input must be a vector')
elseif x.order ~= 1
    error('first input must be a vector')
end

sm = x.size;
st = size(x.basis, 1)*ones(1, x.order);
c1 = reshape(x.components, [sm st]);

if n == 1
    s = 0;
    for i = 1:size(c1,3)
        s = s + abs(c1(:,:,i));
    end
elseif n == 2
    s = 0;
    for i = 1:size(c1,3)
        s = s + c1(:,:,i).^2;
    end
    s = sqrt(s);
else
    s = abs(c1(:,:,1));
    for i = 2:size(c1,3)
        st = abs(c1(:,:,i));
        for j = 1:size(st,1)
            for k = 1:size(st,2)
                if s(j,k) < st(j,k)
                    s(j,k) = st(j,k);
                end
            end
        end
    end
end

end