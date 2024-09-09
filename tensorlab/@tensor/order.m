function o = order(t)
% ORDER  Tensor order.
%  o = ORDER ( X )
%    Returns the order of X
%    o = 1: X is a vector
%    o = 2: X is a second-order tensor
%    o = 3: X is a third-order tensor
%    etc.
% 

if ~isa(t,'tensor')
    o = 0;
else
    o = t.order;
end
