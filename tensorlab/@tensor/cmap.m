function [cm S] = cmap(N,map)

%CMAP define colors for quiver
%     CMAP(X,MAP) maps the lengths of X on the selected colormap MAP.

Nc = N.components;
Nm = N.size;
No = N.order;
Nb = size(N.basis, 1);

S = zeros(1,size(N,1)*No);

X.type = '()';
jj = 1;
for ii = 1:size(N,1)
    X.subs = {[ii]};
    Ni = subsref(N,X);
    
    if No == 1
        
        S(ii) = norm(Ni);
        
    elseif No == 2

        e = eig(Ni);
        S(jj:jj+Nb-1) = e';
        jj = jj + Nb;
        
    end
end

ul = max(S);
ll = min(S);

if ul == ll
    ul = ul + 1;
end

cm = {};
eval(['c = colormap(',map,');']);
for i = 1:length(S)
    cm = {cm{:},interp1(linspace(ll,ul,length(c)),c,S(i),'pchip')};
end   

caxis([ll ul])

end