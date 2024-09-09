function [X,conn] = femgrid ( dimen , numel , e , varargin )
% 
% FEMGRID  Make a square mesh.
% 
%    [X,CONN] = FEMGRID ( DIMEN , NUMEL , E ) creates a square mesh in
%    either 2-D or 3-D depending on the basis E, where E is a column of
%    basis vectors, attainable by for instance
%     E = cartesianbasis2d('ex','ey')
%    for a 2-D Cartesian basis.
% 
%    DIMEN specifies the dimension of the square mesh. I.e. of a 2D basis
%      DIMEN = [ width height ]
%    for a 3D basis
%      DIMEN = [ width height depth ]
%    These are the dimensions in x-, y-, and z-direction respectively.
%
%    NUMEL specifies the number of elements in each direction. I.e. for
%    a 2D basis we have
%      NUMEL = [ el_width el_height ]
%    for a 3D basis we have
%      NUMEL = [ el_width el_height el_depth ]
%    These are the number of in x-, y-, and z-direction respectively.
%    
%    The function returns X which is a column of nodal position vectors
%    and CONN which is the connectivity matrix that contains the mesh
%    topology, i.e. for each element the node numbers that belong to
%    that element. Herein, the relative position within the element
%    is determined by the local numbering. For example the following
%    elements in 2-D:
%
%    Linear triangular element:
%           3   
%          / \
%         /   \ 
%        /     \
%       1-------2
%    then CONN = [1 2 3].
%
%    Quadratic triangular element:
%           3   
%          / \
%         6   5 
%        /     \
%       1---4---2
%    then CONN = [1 2 3 4 5 6].    
%
%    Linear quadrilateral element:
%       4-------3
%       |       |
%       |       |
%       |       |
%       1-------2
%    then CONN = [1 2 3 4].
%
%    Quadratic quadrilaterals element, Serendipity:
%       4---7---3
%       |       |
%       8       6
%       |       |
%       1---5---2
%    then CONN = [1 2 3 4 5 6 7 8].
%
%    Quadratic quadrilaterals element, Lagrange:
%       4---7---3
%       |       |
%       8   9   6
%       |       |
%       1---5---2
%    then CONN = [1 2 3 4 5 6 7 8 9].
% 
%    Next we consider an example of a "real" mesh that consists of linear 
%    quadrilateral elements:
% 
%          elements                           nodes
%    ---------------------         21----22----23----24----25      
%    | 13 | 14 | 15 | 16 |          |     |     |     |     | 
%    ---------------------         16----17----18----19----20
%    |  9 | 10 | 11 | 12 |          |     |     |     |     | 
%    ---------------------         11----12----13----14----15 
%    |  5 |  6 |  7 |  8 |          |     |     |     |     | 
%    ---------------------          6---- 7---- 8---- 9----10
%    |  1 |  2 |  3 |  4 |          |     |     |     |     | 
%    ---------------------          1---- 2---- 3---- 4---- 5
% 
%    This yields the following connectivity matrix:
%    CONN =
% 
%         1     2     7     6
%         2     3     8     7
%         3     4     9     8
%         4     5    10     9
%         6     7    12    11
%         7     8    13    12
%         8     9    14    13
%         9    10    15    14
%        11    12    17    16
%        12    13    18    17
%        13    14    19    18
%        14    15    20    19
%        16    17    22    21
%        17    18    23    22
%        18    19    24    23
%        19    20    25    24
% 
%    Note: the definition of the local numbering that is used here
%    corresponds to that of TENSOR/FEMPLOT. Also, this definition
%    coincides with the definition of the lecture notes of the course
%    4A700 by Peerlings.
% 
%    [X,CONN] = FEMGRID ( ... , 'PropertyName' , property ) allows the
%    choice of different element types. We can choose different element
%    families and within these families different types:
%
%    ELEMENTSHAPE
%      { 'QUAD' } | 'TRIA' herein QUAD corresponds to a hexagonal in
%      3-D or a rectangular element in 2-D. TRIA on the other hand in
%      a tetrahedral in 3-D or a triangular element in 2-D.
%
%      
%    ELEMENTTYPE
%      { 'LINEAR' } | 'SERENDIPITY' | 'LAGRANGE' | 'QUADRATIC'
%      this property determines the number of nodes within each
%      element. Depending on the element family, different element-types
%      are defined. Consider the following overview:
%
%      --------------------------------------------------
%      | FAMILY | TYPE        | #NODES 2-D | #NODES 3-D |
%      --------------------------------------------------
%      | QUAD   | LINEAR      |      4     |      8     |
%      | QUAD   | SERENDIPITY |      8     |     20     |
%      | QUAD   | LAGRANGE    |      9     |     27     |
%      | TRIA   | LINEAR      |      3     |      4     |
%      | TRIA   | QUADRATIC   |      6     |     10     |
%      --------------------------------------------------
% 
%    SEE ALSO: TENSOR/FEMPLOT

% -------------------------------------------------------------------------
% check/convert required input
% -------------------------------------------------------------------------
if length(dimen) == 2
    if length(numel) == 2
        if size(e.basis,1) == 2
            width  = dimen(1);
            height = dimen(2);
            elw    = numel(1);
            elh    = numel(2);
            dim    = 2;
        else, error('Basis dimension does not correspond to the number of elements')
        end
    else, error('Number of elements does not correspond to the dimensions')
    end
elseif length(dimen) == 3
    if length(numel) == 3
        if size(e.basis,1) == 3
            width  = dimen(1);
            height = dimen(2);
            depth  = dimen(3);
            elw    = numel(1);
            elh    = numel(2);
            eld    = numel(3);
            dim    = 3;
        else, error('Basis dimension does not correspond to the number of elements')
        end
    else, error('Number of elements does not correspond to the dimensions')
    end
end

% check the basis
if ~checkbasis(e) || ( size(e.basis,1) ~= 2 && size(e.basis,1) ~= 3 )
   error('Incorrect basis')
end

% -------------------------------------------------------------------------
% read/check optional input
% -------------------------------------------------------------------------

% set defaults
fam = 'QUAD';
typ = 'LINEAR';

% read input
nar=size(varargin,2); i=0;
while i<nar, i=i+1;
	if ~ischar(varargin{i}), error('Property must be a string'); 
	elseif strcmpi(varargin{i},'ElementShape'), i=i+1; fam=varargin{i};
	elseif strcmpi(varargin{i},'ElementType'),   i=i+1; typ=varargin{i};
	else, error('Unknown Property'); 
	end
end

% check settings/convert elementype
% QUAD: 'LINEAR'      > 1
% QUAD: 'SERENDIPITY' > 2
% QUAD: 'LAGRANGE'    > 3
% TRIA: 'LINEAR'      > 1
% TRIA: 'QUADRATIC'   > 3
if ~ischar(fam), error('ElementFamily must be a string'); end
if ~ischar(typ), error('ElementType must be a string');   end
if strcmpi(fam,'QUAD')
	if     strcmpi(typ,'LINEAR'),      typ=1;
	elseif strcmpi(typ,'SERENDIPITY'), typ=2;
	elseif strcmpi(typ,'LAGRANGE'),    typ=3;
	else,  error('Undefined ElementType for this ElementFamily');
	end
elseif strcmpi(fam,'TRIA')
	if     strcmpi(typ,'LINEAR'),      typ=1;
	elseif strcmpi(typ,'QUADRATIC'),   typ=3;
	else,  error('Undefined ElementType for this ElementFamily');
	end
else, error('Undefined ElementFamily'); 
end

% -------------------------------------------------------------------------
% construct the elementtype. Always begin with QUAD elements which can later 
% be converted to TRIA elements if requested
% -------------------------------------------------------------------------

% QUAD:LINEAR / TRIA:LINEAR
if typ==1 

    if dim == 2

        % total number of elements and nodes
        elem  = elw * elh;
        nodes = (elw+1) * (elh+1); 
        % connectivity matrix
        conn = zeros(elem,4);
        % build the connectivity
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
        % nodal positions
        x = zeros(nodes,dim);
		% construct the nodal positions
        nn = 0;
        for hh = 0:elh
            for ww = 0:elw
                nn = nn + 1;
                x(nn,:) = [ width/elw*ww  height/elh*hh ];
            end
        end
     
    else
    
        % total number of elements and nodes
        elem   = elw * elh * eld;
        nodesd = (elw+1) * (elh+1); 
        nodes  = (elw+1) * (elh+1) * (eld+1);    
        % connectivity matrix
        conn = zeros(elem,8);
        % build the connectivity matrix
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
        % nodal positions
        x = zeros(nodes,dim);
        % construct the nodal positions
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
    
% QUAD:LAGRANGE / TRIA:QUADRATIC
elseif typ==3
    
    if dim == 2
        
        % total number of elements and nodes
        elem  = elw * elh;
        nodes = (2*elw+1) * (2*elh+1) ; 
        % connectivity matrix
        conn = zeros(elem,9);
        % build the connectivity matrix
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
        % nodal positions
        x = zeros(nodes,dim);
        % construct the nodal positions
        nn = 0;
        for hh = 0:2*elh
            for ww = 0:2*elw
                nn = nn + 1;
                x(nn,:) = [ width/(2*elw)*ww  height/(2*elh)*hh ];
            end
        end
        
    else
    
        % total number of elements and nodes
        elem   = elw * elh * eld;
        nodesd = (2*elw+1) * (2*elh+1); 
        nodes  = (2*elw+1) * (2*elh+1) * (2*eld+1); 
        % connectivity matrix
        conn = zeros(elem,27);
        % build the connectivity matrix
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
        % nodal positions
        x = zeros(nodes,dim);
        % construct the nodal positions
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
    
% QUAD:SERENDIPITY
elseif typ==2
    
    if dim == 2
        
        % total number of elements and nodes
        elem  = elw * elh;
        nodes = (2*elw+1) * (2*elh+1) - elem; 
        % connectivity matrix
        conn = zeros(elem,8);
        % build the connectivity matrix
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
        % nodal positions
        x = zeros(nodes,dim);
        % construct the nodal positions
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
    
        % total number of elements and nodes
        elem   = elw * elh * eld;
        nodesdf = (2*elw+1) * (2*elh+1) - elw * elh; 
        nodesdp = (elw+1) * (elh+1); 
        nodes  = (nodesdf + nodesdp) * eld + nodesdf; 
        % connectivity matrix
        conn = zeros(elem,20);
        % build the connectivity matrix
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
        % nodal positions
        x = zeros(nodes,dim);
        % construct the nodal positions
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

% set the nodal positions
X = zeros(1,e,size(x,1),1);
X.components = reshape(x,size(x,1)*dim,1);

% -------------------------------------------------------------------------
% change the connectivity for element of the TRIA family
% -------------------------------------------------------------------------

if strcmpi(fam,'TRIA')
	if typ==1 && dim==2
		conn=[ conn(:,[1 2 4]) ; 
               conn(:,[2 3 4]) ]; 
	elseif typ==1 && dim==3
		conn=[ conn(:,[1 2 4 5]) ;
               conn(:,[2 4 6 5]) ;
               conn(:,[5 6 8 4]) ;
               conn(:,[2 3 4 6]) ;
               conn(:,[3 7 6 8]) ;
               conn(:,[4 6 8 3]) ];
	elseif typ==3 && dim==2
		conn=[ conn(:,[1 2 4 5 9 8]) ; 
               conn(:,[2 3 4 6 7 9]) ]; 
	elseif typ==3 && dim==3
		conn=[ conn(:,[1 2 4 5  9 25 12 17 21 24]) ; 
               conn(:,[2 4 6 5 25 27 18 21 24 13]) ;
               conn(:,[5 6 8 4 13 26 16 24 27 20]) ;
               conn(:,[2 3 4 6 10 11 25 18 22 27]) ;
               conn(:,[3 7 6 8 19 14 22 23 15 26]) ;
               conn(:,[4 6 8 3 27 26 20 11 22 23]) ];
	end
end

end


