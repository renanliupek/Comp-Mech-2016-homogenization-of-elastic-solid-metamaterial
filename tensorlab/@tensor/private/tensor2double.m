function d = tensor2double(t)

if t.order ~= 0
    error('Tensorial quantity cannot be changed to double.')
else
    d = reshape(t.components, t.size);
end
