function plotmesh(x,conn,varargin)
% PLOTMESH  Make a plot of a mesh. 
%
%     PLOTMESH(x,conn) plots a mesh. The function plot lines between
%     the nodes that are indicated in conn. The coordinates of the nodes are
%     given in x, which are vectors. 
% 
%     Note: default the plotmesh function renumbers the connectivity matrix,
%     such that for a standard element
%        4---3
%        |   |
%        1---2
%     the conn matrix for that element is [1 2 3 4].
% 
%     Using the MESHGRID([1 1],[1 1],e) function the connectivity matrix is 
%     given as [1 2 4 3].
%     Consult <a href = "matlab:help tensor/meshgrid">help tensor/meshgrid</a> for more information.
% 
%     PLOTMESH(conn,coord,'Property',propertyvalue) is the same as previous but
%     now additional options can be customly selected. The options are
% 
%     'Renumber',
%      {'on'} | 'off'
%       to have the conn matrix renumbered
%       default: if not specified the Renumber option is on
% 
%     'ShowNodes',
%      {'on'} | 'off'
%       to show the node numbers
%       default: if not specified the ShowNodes option is off
% 
%     'PositionNodes',
%      {'auto'},'off',vector
%       the node numbers are plotted at the nodal position. This option allows
%       you to shift the positions at which they are plotted slightly, to
%       improve readability. Default, this function shifts these position
%       slightly; but one can also specify a vector (which must be 2D for a 2D
%       plot and 3D for a 3D plot.
% 
%     'ShowElements',
%      {'on'} | 'off'
%       to show the element numbers (in a black box)
%       default: if not specified the ShowElements option is off
% 


warning('tensor:new','Please use femplot,\nthis function will be removed in newer functions of tensorlab');

% Set defaults
renumber     = 1;
shownodes    = 0;
inputerror   = 0;
Lagrange     = 0;
posnodes     = 1;
showelements = 0;
ecolor       = 'blue';
fcolor       = 'none';
falpha       = 1;

if size(varargin,2) > 12 || mod(size(varargin,2),2) == 1 
    error('Incorrect input. Check input, or consult help plotmesh');
    inputerror = 1;
else
    for i = 1:2:size(varargin,2)
        if strcmpi(varargin(i),'Renumber')
            if strcmp(varargin(i+1),'off')
                renumber = 0;
            elseif ~strcmp(varargin(i+1),'on')
                warning('Tensor:IncorrectInput',...
               'Renumber should be set to "on" or "off" \n Assuming Renumber "on"')
            end
        elseif strcmpi(varargin(i),'ShowNodes')
            if strcmp(varargin(i+1),'on')
                shownodes = 1;
            elseif ~strcmp(varargin(i+1),'off')
                warning('Tensor:IncorrectInput',...
                'ShowNodes should be set to "on" or "off" \n Assuming ShowNodes "off"')
            end
        elseif strcmpi(varargin(i),'PositionNodes')
            if strcmp(varargin(i+1),'off')
                posnodes = 0;
            elseif isa(varargin{i+1},'tensor')
                posnodes = varargin{i+1};
            elseif ~strcmp(varargin(i+1),'auto')
                warning('Tensor:IncorrectInput',...
                'PositionNodes should be set to "auto", "off", or a vector \n PositionNodes ShowNodes "off"')
            end
        elseif strcmpi(varargin(i),'ShowElements')
            if strcmp(varargin(i+1),'on')
                showelements = 1;
            elseif ~strcmp(varargin(i+1),'off')
                warning('Tensor:IncorrectInput',...
                'ShowElements should be set to "on" or "off" \n Assuming ShowNodes "off"')
            end
        elseif strcmpi(varargin(i),'Color') || strcmpi(varargin(i),'EdgeColor')
            ecolor = varargin{i+1};
        elseif strcmpi(varargin(i),'FaceColor')
            fcolor = varargin{i+1};
        elseif strcmpi(varargin(i),'FaceAlpha')
            falpha = varargin{i+1};
        else
            error('Incorrect input. Check input, or consult help plotmesh');
            inputerror = 1;
        end    
    end
end

if ~inputerror
    
dim = size(x.basis,1);

% Renumber the conn matrix (default renumber = 1)
if renumber
    if dim == 2
        if size(conn,2) == 3 
            connplot = conn;
            connex = [];
        elseif size(conn,2) == 6
            connplot = conn(:,1:3);
            connex = [];
        elseif size(conn,2) == 4 % 2D linear elements
            connplot = conn;
            connex = [];
        elseif size(conn,2) == 8 || size(conn,2) == 9 % 2D Serendipity or Lagrange
            connplot = conn;
            connplot(:,1:9) = connplot(:,[1 5 2 6 3 7 4 8 1]);
            
             % The midpoints are not used to construct the planes (only
             % Lagrange)
            if size(conn,2) == 8;
                connex = [];
            else
                Lagrange = 1;
                connex = conn(:,9);
            end
        else
            warning('Tensor:connectivity','Unknown dimension connectivity matrix')
        end
        
    elseif dim == 3
        if size(conn,2) == 4 % 3D triangular elements
            connplot = [ conn(:,[1 2 3]) 
                         conn(:,[1 2 4]) 
                         conn(:,[1 3 4]) 
                         conn(:,[2 3 4]) ];
             connex = [];
        elseif size(conn,2) == 8 % 3D linear elements
            connplot = zeros(6*size(conn,1),5);
            jj = 1;
            for ii = 1:size(conn,1)
                connplot(jj,:) = conn(ii,[ 1 2 3 4 1 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[ 1 2 6 5 1 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[ 2 6 7 3 2 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[ 3 7 8 4 3 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[ 4 8 5 1 4 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[ 5 6 7 8 5 ]); jj = jj + 1;
            end
            connex = [];
        elseif size(conn,2) == 20 || size(conn,2) == 27 % 3D Serendipity or Lagrange
            connplot = zeros(6*size(conn,1),9);
            jj = 1;
            for ii = 1:size(conn,1)
                connplot(jj,:) = conn(ii,[  1  9  2 10  3 11  4 12  1 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[  5 13  6 14  7 15  8 16  5 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[  1  9  2 18  6 13  5 17  1 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[  2 10  3 19  7 14  6 18  2 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[  3 11  4 20  8 15  7 19  3 ]); jj = jj + 1;
                connplot(jj,:) = conn(ii,[  4 20  8 16  5 17  1 12  4 ]); jj = jj + 1;
            end
            
            % The midpoints are not used to construct the planes (only
             % Lagrange)
            if size(conn,2) == 27
                Lagrange = 1;
                connex = zeros(size(conn,1),7);
                for ii = 1:size(conn,1)
                    connex(ii,:)   = conn(ii,[ 25 26 27 21 22 23 24 ]);
                end
            else
                connex = [];
            end
            
        else
            warning('Tensor:connectivity','Unknown dimension connectivity matrix')
        end
    else
        error('Incorrect dimension')
    end
end

% Sweep the connplot for any double entries. If this would not be done, all
% these planes would be drawn twice.
connplot = unique(connplot,'rows');
% remove = [];
% for ii = 1:size(connplot,1)
%     for jj = 1:size(connplot,1)
%         if ii > jj
%             bool = connplot(ii,:) == connplot(jj,:);
%             i = 1;
%             while i < size(connplot,2) && bool(i) == 1
%                 i = i + 1;
%             end
%             if i < size(connplot,2)
%                 bool = 0;
%             else
%                 bool = bool(i);
%             end
%             
%             if bool
%                 remove = [remove ii];
%             end
%         end
%     end
% end
% 
% % Keep only the non-double parts
% connplot = connplot(setdiff(1:size(connplot,1),remove),:);

% Sweep the 'extra' points that need to be plotted for the Lagrange
% elements, for any double entries.

% remove = [];
connex = reshape(connex,size(connex,1)*size(connex,2),1);
connex = unique(connex);
% for ii = 1:size(connex,1)
%     for jj = 1:size(connex,1)
%         if ii > jj
%             if connex(ii) == connex(jj)
%                 remove = [remove ii];
%             end
%         end
%     end
% end
% 
% % Keep only the non-double parts
% connex = connex(setdiff(1:size(connex,1),remove),:);    

% Make the plot

if ~shownodes
    
    patch('Faces',connplot,'Vertices',x,'EdgeColor',ecolor,'FaceColor','none',...
          'FaceColor',fcolor,'FaceAlpha',falpha)
    
else
    
    patch('Faces',connplot,'Vertices',x,'EdgeColor',ecolor,'FaceColor','none',...
          'FaceColor',fcolor,'FaceAlpha',falpha,...
          'Marker','o','MarkerEdgeColor',ecolor,'MarkerFaceColor',ecolor)
                
end

% For a 3D plot, reset the view
if dim == 3
    view(45,45);
end

% Reset the axis of the figure
X = reshape(x.components, size(x.components,1)/size(x.basis,1), size(x.basis,1) );

ll = min(X);
ul = max(X);

fac = []; % Shape factor
ax = reshape([ll;ul],1,2*dim);
for ii = 1:2:2*dim
    sh = abs(ax(ii)-ax(ii+1));
    ax(ii)   = ax(ii)   - .1*sh;
    ax(ii+1) = ax(ii+1) + .1*sh;
    fac = [fac sh];
end

ac = axis;

if size(ax,2) ~= size(ac,2)
    error('Close the figure first, before changing dimension\n');
end

ll = min([ax;ac]);
ul = max([ax;ac]);

for ii = 1:2:2*dim
    ac(ii)   = ll(ii);
    ac(ii+1) = ul(ii+1);
end

axis(ac)

% Plot the node numbers

if shownodes
    
    S.type = '()';
    
    hold on

    for ii = 1:size(x,1)
        
        node = sprintf('%d',ii);
        
        S.subs = {[ii]};
        xcv = subsref(x,S);
        
        if ~isa(posnodes,'tensor')
            if posnodes
                xcv.components = xcv.components + .01*fac';
            end
        elseif isa(posnodes,'tensor')
            xcv.components = xcv.components + posnodes.components;
        end
        
        text(xcv,node,'Color',ecolor);
        
    end
    
    for ii = 1:size(connex,1)
                
        S.subs = {[connex(ii)]};
        xii = subsref(x,S);
        xii = [xii.components];

        if size(xii,1) == 2
            plot(xii(1),xii(2),'LineStyle','none','Color',ecolor,...
                'Marker','o','MarkerEdgeColor',ecolor,'MarkerFaceColor',ecolor)
        elseif size(xii,1) == 3
            plot3(xii(1),xii(2),xii(3),'LineStyle','none','Color',ecolor,...
                'Marker','o','MarkerEdgeColor',ecolor,'MarkerFaceColor',ecolor)
        end
                    
    end
            
    hold off
    
end

% Plot the element numbers

if showelements
    
    hold on
            
    for i = 1:size(conn,1)

        elem = sprintf('%d',i);

        S.type = '()';
        S.subs = {[conn(i,:)]};
        
        xiv = subsref(x,S);
        xiv = mean(xiv);

        if Lagrange
            xiv.components = xiv.components - .05*fac';
        end
        
        text(xiv,elem,'Color',ecolor,'EdgeColor',ecolor)

    end
    
    hold off
    
end


end


end
