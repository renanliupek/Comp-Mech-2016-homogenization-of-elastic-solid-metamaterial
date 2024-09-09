function [ I ] = identity ( order ,e )

% IDENTITY  Identity tensor. 
% 
%   I = IDENTITY ( order , e )
%   Returns the identity tensor or order (order). To construct a tensor, a
%   basis is always needed. (e) is therefore a column of basis vectors.
% 
%   Notice that the identity tensor can by definition only be of an even
%   order.

if ~isa( order , 'float' )
    error('The order must be numeric')
elseif mod(order,2)
    error('Only even order identity tensors exist')
elseif order > 16
    error('The identity tensor is not defined for an order higher than 16')
elseif order < 2
    error('The identity tensor is only defined for order 2 and higher')
elseif ~isa( e , 'tensor' )
    error('The basis must be a tensor object')
elseif ~checkbasis( e )
    error('Incorrect basis input')
end

dim = size(e.basis,1);

if order == 2
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
        c(:,:,i,i) = 1;
    end
elseif order == 4
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
        for j = 1:dim
                c(:,:,i,j,j,i) = 1;
        end
    end
elseif order == 6
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       c(:,:,i,j,k,k,j,i) = 1;
      end
     end
    end
elseif order == 8
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       for l = 1:dim
        c(:,:,i,j,k,l,l,k,j,i) = 1;
       end
      end
     end
    end
elseif order == 10
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       for l = 1:dim
        for m = 1:dim
         c(:,:,i,j,k,l,m,m,l,k,j,i) = 1;
        end
       end
      end
     end
    end
elseif order == 12
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       for l = 1:dim
        for m = 1:dim
         for n = 1:dim
          c(:,:,i,j,k,l,m,n,n,m,l,k,j,i) = 1;
         end
        end
       end
      end
     end
    end
elseif order == 14
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       for l = 1:dim
        for m = 1:dim
         for n = 1:dim
          for o = 1:dim
           c(:,:,i,j,k,l,m,n,o,o,n,m,l,k,j,i) = 1;
          end
         end
        end
       end
      end
     end
    end
elseif order == 16
    c = zeros([1 1 dim*ones(1,order)]);
    for i = 1:dim
     for j = 1:dim
      for k = 1:dim
       for l = 1:dim
        for m = 1:dim
         for n = 1:dim
          for o = 1:dim
           for p = 1:dim
            c(:,:,i,j,k,l,m,n,o,p,p,o,n,m,l,k,j,i) = 1;
           end
          end
         end
        end
       end
      end
     end
    end
end

I.basis = e.basis;
I.order = order;
I.size = [1 1];
I.components = c(:);

I = tensor(I);

end
