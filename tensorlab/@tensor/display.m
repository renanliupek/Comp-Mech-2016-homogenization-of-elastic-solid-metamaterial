function display(t)

% echo variable name
fprintf('\n%s =\n\n', inputname(1))

% empty matrix
if isempty(t.components)

    fprintf('     []\n')

% scalar matrix
elseif t.order == 0
    
    disp(reshape(t.components, t.size));

% tensors and matrices of tensors
else
    
    % pre-factor for small and large values
    m = max(abs(t.components));
    e = 0;
    if (m >= 1e-16) && ( (m < 1e-3) || (m > 1e4) )
        e = fix(log10(m));
        fprintf('   %7.1e *\n\n', 10^e)
        t.components = t.components * 10^-e;
    end
    
    % loop over matrix elements
    mm = ones(1, length(t.size));
    mm(1) = 0;
    for m = 1:prod(t.size)
        
        % construct matrix indices
        mm(1) = mm(1) + 1;
        for n = 1:length(t.size)
            if mm(n) > t.size(n)
                mm(n) = 1;
                mm(n+1) = mm(n+1) + 1;
            end
        end
        
        % print matrix indices except for single column
        if max(t.size(2:end)) > 1
            fprintf('   (');
            for n = 1:length(t.size)
                fprintf(sprintf('%%%dd', length(num2str(t.size(n)))), mm(n))
                if n < length(t.size)
                    fprintf(',')
                else
                    fprintf(')  ')
                end
            end
        end
        
        % loop over dyadic terms of matrix element
        ii = ones(1, t.order);
        ii(end) = 0;
        for i = 1:size(t.basis, 1)^t.order
            
            % construct dyads
            ii(end) = ii(end) + 1;
            for j = t.order:-1:1
                if ii(j) > size(t.basis, 1)
                    ii(j) = 1;
                    ii(j-1) = ii(j-1) + 1;
                end
            end
            
            % retrieve coefficient
            j = m + prod(t.size) * sum((ii-1).*size(t.basis, 1).^(0:t.order-1));
            
            % print term
            if i == 1
                if t.components(j) < 0
                    fprintf('   -')
                else
                    fprintf('    ')
                end
            else
                if t.components(j) < 0
                    fprintf(' - ')
                else
                    fprintf(' + ')
                end
            end
            fprintf('%6.4f', abs(t.components(j)))
            for j = 1:t.order
                fprintf('*%s', deblank(t.basis(ii(j), :)))
            end

        end

        % finish line
        fprintf('\n')
            
    end
    
end

% finish display
fprintf('\n')
