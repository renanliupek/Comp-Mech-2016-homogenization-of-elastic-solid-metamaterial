function q=mesh2d_lin2qua_uc(l)

% INPUT: struct with comsol linear mesh
% OUTPUT: struct with comsol quadratic mesh

% relies on tensorlab



% INPUT: cell/ container.map with linear mesh data
% Notation from TensorLab
% Linear quadrilateral element:
%         4-------3
%         |       |
%         |       |
%         |       |
%         1-------2
%      then CONN = [1 2 3 4].

% OUTPUT: cell/ container.map with quadratic mesh data
% Notation from TensorLab
% Quadratic serendipity element:
%         4----7----3
%         |         |
%         8         6
%         |         |
%         1----5----2
%      then CONN = [1 2 3 4 5 6 7 8].



% %% ADD TENSORLAB & SHARED FUNCTIONS & MESHES TO PATH
% path(path,['C:\GitLab\Tools\TensorLab20\TensorLab20']);
% path(path,["C:\GitLab\Tools\TensorLab20\TensorLab20"+...
%             "\dev\shared_functions"]);
% path(path,["C:\GitLab\Tools\TensorLab20\TensorLab20"+...
%             "\dev\meshes"]);

%% DEFINE TENSORLAB BASIS
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(b{1}, b{2});
e1 = ee(1);
e2 = ee(2);

%% MANUALLY IMPORT COMSOL MESH AS STRUCT WITH CELL ELEMENTS
%--------------------------------------------------------------------------
% open txt file comsol_mesh_solid_rect.mphtxt and copy/paste the node 
% position and connectivity matrix.
% -----------------------------x(--------------------------------------------
%  Comments on 'conn' and 'tag' matrices for elements and edges      STATUS
% -------------------------------------------------------------------------
% - swap column 3 and 4 of conn ( quadrilateral element);                OK
% - no need for swaping when triangular element;                         OK
% - add 1 to all nodes of 'conn' (maybe 'tag' too). In comsol it starts
% in 0, whereas in matlab, it starts in 1;                               OK
% - attention to units! Use SI.                                          OK
% - check if the mesh has at least 6 nodes per wavelength          ON-GOING
%
%
%
% -------------------------------------------------------------------------
%  Comments on edge/element tags in comsol mesh file                                
% -------------------------------------------------------------------------
% - list of fluid-solid edge tag: 4,6
% This means that if an edge has tag 4 or 6, than the edge is at the fluid-
% solid interface.
% - list of element type (fluid or solid) after 'conn' matrix
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

% DESCRIPTION OF STRUCT MESH 2 X 2 
% l.elem =
%            'conn'	            'tag'
%         lm x 2 double	     lm x 1 double
%
% l.edge =
%            'conn'	            'tag'
%         lm x 2 double	     lm x 1 double
%
% l.x =
%         'coordinate'	       'tensor vec'
%         ln x 2 double	     ln x 1 double
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOAD : desactivated when script becomes a function
% % load manually imported comsol mesh
% load('es_fluidsolid_rect_MESH.mat');
% 
% % linear element mesh
% l=mesh;
% clear mesh;
%--------------------------------------------------------------------------




% mesh parameters, linear elements
m   = size(l.elem{2,1}, 1); % number of volume elements in the mesh
n   = size(l.x{2,2}   , 1); % number of nodes in the mesh 
m_e = size (l.edge{2,1},1); % number of edge elements in the mesh


%% INITIALIZE QUADRATIC MESH

% define the cell content with strings
q.elem = {'conn', 'tag'};
q.edge = {'conn', 'tag'};
q.x    = {'coordinate', 'tensor vec'};


% VOLUME ELEMENTS q.elem
 % conn ( add 4 columns)
 q.elem{2,1}         = zeros([m 8]);
 q.elem{2,1}(:,1:4)  = l.elem{2,1};
 % tag (same as linear mesh)
 q.elem{2,2} = l.elem{2,2};

% EDGE ELEMENTS q.edge
  if isempty(l.edge{2,1}); error("incorrect mesh format");end
  % # of elements 
  q.edge{2,1}              =  zeros(m_e,3);
  q.edge{2,1}(:,1)         =  l.edge{2,1}(:,1);
  q.edge{2,1}(:,3)         =  l.edge{2,1}(:,2);
  % tag  (same as linaer)
  q.edge{2,2}              = l.edge{2,2}  ;
  

% POSITION
 % declare 'tensor vec'
%  q.x{2,2} = zeros(1, ee,  n,1); % 1-order tensor
 % assign
 q.x{2,2} = l.x{2,2}; 


%% FOR LOOP: BUILD CONN OF ALL EDGE ELEMENTS
% apparently, comsol conn matrix for edge elements dont have all the edges
% in it. Using 'for loop' on the volume elements, stacking the edges in 
% 'conn_ed_all', and after using 'unique' allows to get rid of duplicates.
conn_ed_all=[];
for e=1:m
    % extract initial element data
    iie = l.elem{2,1}(e,:) ;    % nodes of the current element e
%     xe  = x(iie) ;     % coordinates of these nodes
    conn_ed_all = [conn_ed_all;
                   iie(1) iie(2);
                   iie(2) iie(3);
                   iie(3) iie(4);
                   iie(4) iie(1)];
end

% eliminate duplicates with reversed node order
conn_ed_all = unique( conn_ed_all, 'rows');
emax=size(conn_ed_all,1);
e=1;
while e<emax
    % swap nodes positions of column e
    conn_ed_all(e,[1 2])=conn_ed_all(e,[2 1]);
    % compute new conn without duplicates and store in aux
    conn_ed_all_aux = unique(conn_ed_all, 'rows');
    % if there was a duplicate edge enter if below
    if size(conn_ed_all_aux,1)<size(conn_ed_all,1)
        conn_ed_all = conn_ed_all_aux;
        emax=size(conn_ed_all,1);
%         e=1;   % this is not optimized choice.
        if e~= 1; e=e-1;end % go back one element to check duplicates

    else
    e=e+1;
    end
end

m_e_all = size(conn_ed_all,1);


%% FOR LOOP: ADD MIDDLE POINT NODES VOLUME ELEMENT

% initial conn matrix of quadratric elements
conn = q.elem{2,1};
% tensor vector
x    = q.x{2,2};
% vx   = q.x{2,1};

% initial number of nodes
node_count = n; 

% conn of edge elements (COMSOL)
conn_ed = q.edge{2,1};

%--------------------------------------------------------------------------
% Before loop, master element
%               edge element
%           1------------------2
%
% After loop, master element
%               edge element
%           1--------3---------2
%
%--------------------------------------------------------------------------
for e = 1:m_e_all % loop over mesh all edge elements

% EDGE ELEMENT
    % extract earliest edge data (still with second column of zeros)
    iie = conn_ed_all(e, :) ;    % nodes of the current element e
    xe  = x(iie) ;        % coordinates of these nodes
    % add a node
    node_count     = node_count + 1;
%     iie(2)         = node_count;
    x(node_count)  = (xe(1) + xe(2))/2;              

   
    % assign middle node to comsol edge conn
    % find element 'ec' in comsol mesh 
    [el1,~] = find(conn_ed==iie(1));
    [el2,~] = find(conn_ed==iie(2));
    ec     = intersect(el1,el2);
    if ~isempty(ec)
    conn_ed(ec,2)=node_count;
    end

% VOLUME ELEMENT
    % the edge 'e' is and edge of 1 volume element ( at the boundary) or it
    % is an edge of 2 volume elements (inner domain)
    % 1-3-2 (edge element)
    [row1,col1]=find(conn==iie(1)); % node 1
    [row2,col2]=find(conn==iie(2)); % node 2
    is_edge_of_elems = intersect(row1,row2);
    for j=1:length(is_edge_of_elems)
        % find which nodes in the volume element this edge corresponds
        %nodes(1)=col1(find(is_edge_of_elems(j)==row1));
        %nodes(2)=col2(find(is_edge_of_elems(j)==row2));
        %use logical array as index
        nodes(1)=col1((is_edge_of_elems(j)==row1));
        nodes(2)=col2((is_edge_of_elems(j)==row2));
        if all(ismember(nodes,[1 2])) % if 1-5-2 edge
        conn(is_edge_of_elems(j),5)=node_count;
        end
        if all(ismember(nodes,[2 3])) % if 2-6-3 edge
        conn(is_edge_of_elems(j),6)=node_count;
        end
        if all(ismember(nodes,[3 4])) % if 3-7-4 edge
        conn(is_edge_of_elems(j),7)=node_count;
        end
        if all(ismember(nodes,[4 1])) % if 4-8-1 edge
        conn(is_edge_of_elems(j),8)=node_count;
        end
    end
end

q.edge{2,1} = conn_ed;
q.elem{2,1} = conn;
q.x{2,2}    = x;
% clear conn;

%% PLOT QUADRATIC MESH
% figure(2)
% clf;
% femplot(q.x{2,2} ,q.elem{2,1},'Color' , 'black', 'Nodes', 'off', ...
% 'NodeNumbers' , 'on', 'ElementNumbers', 'on');
% title('quadratic element mesh')

%% RENAMING NODES
% It is interesting to have the node number relatively close within each
% element because the mass/stiff matrices have nonzero elements closer to
% the diagonal. That's called a 'narrow' band for a sparse matrices.


%--------------------------------------------------------------------------
% unit cell size in x direction
ax = 2*max(dot(q.x{2,2},e1));

% unit cell size in y direction
ay = 2*max(dot(q.x{2,2},e2));
%--------------------------------------------------------------------------

% get the index linking list, assuming the center of the unit cell is at
% the point (0,0). The node 1 will be at position (-a/2,-a/2).
[~,linklist] = sort(norm(x+ (ax/2*e1+ay/2*e2)*ones(size(x,1),1) ));

% declare new matrices to be modified
newconn    = conn;
newconn_ed = conn_ed;

% modify new matrices according to linking list
for i=1:length(linklist)
    % substitute node linklist(i) by node i in matrices conn and conn_ed
    newconn   (conn    ==linklist(i))=i;
    newconn_ed(conn_ed ==linklist(i))=i;
end

% assign to mesh, note tags remain the same regardless of node change
q.elem{2,1}   = newconn;
q.edge{2,1}   = newconn_ed;
q.x{2,2}      = q.x{2,2}(linklist);
q.x{2,1}(:,1) = dot(q.x{2,2},e1);
q.x{2,1}(:,2) = dot(q.x{2,2},e2);


%% PLOT FINAL QUADRATIC MESH
% figure(3)
% clf;
% femplot(q.x{2,2} ,q.elem{2,1},'Color' , 'black', 'Nodes', 'off', ...
% 'NodeNumbers' , 'on', 'ElementNumbers', 'on');
% title('quadratic element mesh (renumbered nodes)')

 %% SAVE QUADRATIC ELEMENT MESH
% 
% save([pwd '\meshes\es_fluidsolid_rect_QUADRATIC_MESH.mat'], "q");
end