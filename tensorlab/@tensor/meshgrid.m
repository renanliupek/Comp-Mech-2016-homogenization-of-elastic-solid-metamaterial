function [X,conn] = meshgrid ( dimen , numel , e , varargin )
% 
% MESHGRID  Make a square mesh.
% 
%     [x,conn] = MESHGRID ( DIMEN , NUMEL , E )
%     creates a square mesh in either 2D or 3D depending on the basis E. Here E
%     is a column of basis vectors, attainable by for instance
%     E = cartesianbasis2d('ex','ey')
% 
%     DIMEN specifies the dimension of the square mesh. I.e. of a 2D basis
%      DIMEN = [ width height ]
%     for a 3D basis
%      DIMEN = [ widht height depth ]
%     These are the dimensions in x-, y-, and z-direction respectively.
% 
%     NUMEL specifies the number of elements in each direction. I.e. for a 2D
%     basis we have
%      NUMEL = [ el_widht el_height ]
%     for a 3D basis we have
%      NUMEL = [ el_widht el_height el_depth ]
%     These are the number of in x-, y-, and z-direction respectively.
% 
%     The function returns x which is a column of nodal position vectors and
%     conn which is the connectivity matrix that contains the mesh topology,
%     i.e. per element the nodes of that element
% 
%     For example, consider the following mesh
% 
%     elements              nodes
%     13 14 15 16           21 22 23 24 25
%      9 10 11 12           16 17 18 19 20
%      5  6  7  8            6  7  8  9 10
%      1  2  3  4            1  2  3  4  5
% 
%     This yield the conn matrix
%     conn =
% 
%          1     2     7     6
%          2     3     8     7
%          3     4     9     8
%          4     5    10     9
%          6     7    12    11
%          7     8    13    12
%          8     9    14    13
%          9    10    15    14
%         11    12    17    16
%         12    13    18    17
%         13    14    19    18
%         14    15    20    19
%         16    17    22    21
%         17    18    23    22
%         18    19    24    23
%         19    20    25    24
% 
%     Note: if you use this conn matrix in plotmesh, the option 'Renumber'
%     should be set to 'on' (default).
% 
%     Note2: the order of entries in the conn matrix is the standard numbering
%     of a single element (see Peerlings).
% 
%     [x,conn] = MESHGRID ( DIMEN , NUMEL , E , 'PropertyName',property)
%     is the same as previous, but now allows to choose different element
%     types. The choices are
%      'ElementType','Linear'      for 4 noded elements in 2D, and 8 in 3D
%      'ElementType','Serendipity' for 8 noded elements in 2D, and 20 in 3D
%      'ElementType','Lagrange'    for 9 noded elements in 2D, and 27 in 3D
% 
%     To plot the mesh, use <a href = "matlab:help tensor/plotmesh">plotmesh</a>

warning('tensor:new','Please use femgrid,\nthis function will be removed in newer functions of tensorlab');

if length(dimen) == 2
    if length(numel) == 2
        if size(e.basis,1) == 2
            width = dimen(1);
            height = dimen(2);
            elw = numel(1);
            elh = numel(2);
            dim = 2;
        else
            error('Basis dimension does not correspond to the number of elements')
        end
    else
        error('Number of elements does not correspond to the dimensions')
    end
elseif length(dimen) == 3
    if length(numel) == 3
        if size(e.basis,1) == 3
            width = dimen(1);
            height = dimen(2);
            depth = dimen(3);
            elw = numel(1);
            elh = numel(2);
            eld = numel(3);
            dim = 3;
        else
            error('Basis dimension does not correspond to the number of elements')
        end
    else
        error('Number of elements does not correspond to the dimensions')
    end
end

Linear = 0;
Lagrange = 0;
Serendipity = 0;
elementtype = 0;

i = 1;
while i < size(varargin,2)
    if strcmp(varargin{i},'ElementType')
        elementtype = 1;
        i = i + 1;
        if strcmp(varargin{i},'Linear')
            Linear = 1;
            i = i + 1;
        elseif strcmp(varargin{i},'Lagrange')
            Lagrange = 1;
            i = i + 1;
        elseif strcmp(varargin{i},'Serendipity')
            Serendipity = 1;
            i = i + 1;
        else
            error('Unknown element type')
        end
        
    else
        
        error('Unknown input');
    end
end

% Check the basis
if ~checkbasis(e) || ( size(e.basis,1) ~= 2 && size(e.basis,1) ~= 3 )
   error('Incorrect basis')
end

% Define a default element (quadrilateral) is non is supplied
if ~elementtype
    Linear = 1;
    warning('Assuming linear element')
end


if Linear
    if dim == 2
        
    % The total number of elements and nodes
    elem  = elw * elh;
    nodes = (elw+1) * (elh+1); 
    
    % Connectivity matrix
    conn = zeros(elem,4);
    
    el = 0;
    for hh = 0:elh-1
        for bb = 0:elw-1
            el = el + 1;
            conn(el,:) = hh*(elw+1)*ones(1,4) + ... 
                         [                  1+bb 2+bb ...
                          (elw+1)*[1 1] + [ 1+bb 2+bb ] ];
        end
    end
    
    conn = conn(:,[1 2 4 3]);
    
    % Node coordinates
    x = zeros(nodes,dim);
    
    nn = 0;
    for hh = 0:elh
        for ww = 0:elw
            nn = nn + 1;
            x(nn,:) = [ width/elw*ww  height/elh*hh ];
        end
    end
        
    else
    
    % The total number of elements and nodes
    elem   = elw * elh * eld;
    nodesd = (elw+1) * (elh+1); 
    nodes  = (elw+1) * (elh+1) * (eld+1); 
    
    % Connectivity matrix
    conn = zeros(elem,8);
    
    el = 0;
    for dd = 0:eld-1
        for hh = 0:elh-1
            for bb = 0:elw-1
                el = el + 1;
                conn(el,:) = hh*(elw+1)*ones(1,8) + dd*nodesd*ones(1,8) + ... 
                             [                      1+bb 2+bb ...
                              (elw+1)*ones(1,2) + [ 1+bb 2+bb ] ...
                               nodesd*ones(1,4) + [ 1+bb 2+bb ...
                              (elw+1)*ones(1,2) + [ 1+bb 2+bb ] ] ];
            end
        end
    end
    
    conn = conn(:,[1 2 4 3 5 6 8 7]);
    
    % Node coordinates
    x = zeros(nodes,dim);
    
    nn = 0;
    for dd = 0:eld
        for hh = 0:elh
            for ww = 0:elw
                nn = nn + 1;
                x(nn,:) = [ width/elw*ww  height/elh*hh  depth/eld*dd ];
            end
        end
    end
    
    end
    
elseif Lagrange
    
    if dim == 2
        
    % The total number of elements and nodes
    elem  = elw * elh;
    nodes = (2*elw+1) * (2*elh+1) ; 
    
    % Connectivity matrix
    conn = zeros(elem,9);
    
    el = 0;
    for hh = 0:elh-1
        for bb = 0:elw-1
            el = el + 1;
            conn(el,:) = hh*(4*elw+2)*ones(1,9) + ... 
                         [                        1+2*bb 2+2*bb 3+2*bb ...
                          (2*elw+1)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ...
                          (4*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ];
        end
    end
    
    conn(:,1:9) = conn(:,[1 3 9 7 2 6 8 4 5]);
    
    % Node coordinates
    x = zeros(nodes,dim);
    
    nn = 0;
    for hh = 0:2*elh
        for ww = 0:2*elw
            nn = nn + 1;
            x(nn,:) = [ width/(2*elw)*ww  height/(2*elh)*hh ];
        end
    end
        
    else
    
    % The total number of elements and nodes
    elem   = elw * elh * eld;
    nodesd = (2*elw+1) * (2*elh+1); 
    nodes  = (2*elw+1) * (2*elh+1) * (2*eld+1); 
    
    % Connectivity matrix
    conn = zeros(elem,27);
    
    el = 0;
    for dd = 0:eld-1
        for hh = 0:elh-1
            for bb = 0:elw-1
                el = el + 1;
                conn(el,:) = hh*(4*elw+2)*ones(1,27) + dd*2*nodesd*ones(1,27) + ... 
                 [                        1+2*bb 2+2*bb 3+2*bb ...
                  (2*elw+1)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ...
                  (4*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ...
                     nodesd*ones(1,9) + [ 1+2*bb 2+2*bb 3+2*bb ...
                  (2*elw+1)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ...
                  (4*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ] ...
                   2*nodesd*ones(1,9) + [ 1+2*bb 2+2*bb 3+2*bb ...
                  (2*elw+1)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ...
                  (4*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ] ];
            end
        end
    end
    
    conn = conn(:,[1 3 9 7 19 21 27 25 2 6 8 4 20 24 26 22 10 12 18 16 11 15 17 13 5 23 14]);
    
    x = zeros(nodes,dim);
    
    nn = 0;
    for dd = 0:2*eld
        for hh = 0:2*elh
            for ww = 0:2*elw
                nn = nn + 1;
                x(nn,:) = [ width/(2*elw)*ww  height/(2*elh)*hh  depth/(2*eld)*dd ];
            end
        end
    end
    
    end
    

elseif Serendipity
    
    if dim == 2
        
    % The total number of elements and nodes
    elem  = elw * elh;
    nodes = (2*elw+1) * (2*elh+1) - elem; 
    
    % Connectivity matrix
    conn = zeros(elem,8);
    
    el = 0;
    for hh = 0:elh-1
        for bb = 0:elw-1
            el = el + 1;
            conn(el,:) = hh*(3*elw+2)*ones(1,8) + ... 
                         [                        1+2*bb 2+2*bb 3+2*bb ...
                          (2*elw+1)*ones(1,2) + [ 1+  bb 2+  bb        ] ...
                          (3*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ];
        end
    end
    
    conn = conn(:,[1 3 8 6 2 5 7 4]);
    
    % Node coordinates
    x = zeros(nodes,dim);
    
    nn = 0;
    for hh = 0:2*elh
        if mod(hh,2) == 0
            for ww = 0:2*elw
                nn = nn + 1;
                x(nn,:) = [ width/(2*elw)*ww  height/(2*elh)*hh ];
            end
        else
            for ww = 0:elw
                nn = nn + 1;
                x(nn,:) = [ width/(elw)*ww    height/(2*elh)*hh ];
            end
        end
    end
        
    else
    
    % The total number of elements and nodes
    elem   = elw * elh * eld;
    nodesdf = (2*elw+1) * (2*elh+1) - elw * elh; 
    nodesdp = (elw+1) * (elh+1); 
    nodes  = (nodesdf + nodesdp) * eld + nodesdf; 
    
    % Connectivity matrix
    conn = zeros(elem,20);
    
    el = 0;
    for dd = 0:eld-1
        for hh = 0:elh-1
            for bb = 0:elw-1
                el = el + 1;
                conn(el,:) = [ hh*(3*elw+2)*ones(1,8) + dd*(nodesdf+nodesdp)*ones(1,8) + ... 
                 [                        1+2*bb 2+2*bb 3+2*bb     ...
                  (2*elw+1)*ones(1,2) + [ 1+  bb 2+  bb        ]   ...
                  (3*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ]   ...
                 hh*(elw+1)*ones(1,4) + dd*(nodesdf+nodesdp)*ones(1,4) + ...
                    nodesdf*ones(1,4) + [ 1+  bb 2+  bb            ...
                  (  elw+1)*ones(1,2) + [ 1+  bb 2+  bb        ] ] ...
               hh*(3*elw+2)*ones(1,8) + dd*(nodesdf+nodesdp)*ones(1,8) + ... 
          (nodesdf+nodesdp)*ones(1,8) + [ 1+2*bb 2+2*bb 3+2*bb     ...
                  (2*elw+1)*ones(1,2) + [ 1+  bb 2+  bb        ]   ...
                  (3*elw+2)*ones(1,3) + [ 1+2*bb 2+2*bb 3+2*bb ] ] ];
            end
        end
    end
    
    conn = conn(:,[1 3 8 6 13 15 20 18 2 5 7 4 14 17 19 16 9 10 12 11]);
    
    x = zeros(nodes,dim);
    
    nn = 0;
    for dd = 0:2*eld
        if mod(dd,2) == 0
            for hh = 0:2*elh
                if mod(hh,2) == 0
                    for ww = 0:2*elw
                        nn = nn + 1;
                        x(nn,:) = [ width/(2*elw)*ww  height/(2*elh)*hh  depth/(2*eld)*dd ];
                    end
                else
                    for ww = 0:elw
                        nn = nn + 1;
                        x(nn,:) = [ width/(elw)*ww    height/(2*elh)*hh  depth/(2*eld)*dd ];
                    end
                end
            end
        else
            for hh = 0:elh
                for ww = 0:elw
                    nn = nn + 1;
                    x(nn,:) = [ width/elw*ww  height/elh*hh  depth/(2*eld)*dd ];
                end
            end
        end
    end
    
    end
              
    
end

X = zeros(1,e,size(x,1),1);
X.components = reshape(x,size(x,1)*dim,1);


end
