function plotdomain(mx,conn, color)
% function plots the mesh 
% figure
hold on

% number of elements
m = size(conn,1);

% plot the mesh
for e =1:m
     
    iie =   conn(e, 1:4); % nodes of the current element e
    xe  =   mx(iie,:);     % coordinates of these nodes
     
     fill(xe(:,1),xe(:,2),color,'edgecolor','none');
     
end

axis equal
axis tight

end
