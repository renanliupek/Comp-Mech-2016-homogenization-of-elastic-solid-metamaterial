function varargout=femmesh(x,t,varargin)
% FEMPLOT  Make a plot of a mesh. 
%
%    FEMPLOT(X,CONN) plots a mesh. The function plots the element in
%    faces. The faces are assigned no color by default. Therefore, element
%    edges appear as lines. The input comprises the nodal positions X (Nx1
%    column of vectors, with N the number of nodes) and the connectivity
%    CONN (MxT matrix of integers, with M the number of elements and T
%    depending on the type of element as well as the dimension). The
%    matrix should be numbered according  the local numbering, which
%    depends on the element type. For 2-D elements, the local numbering
%    is included below. For the 3-D equivalents the literature can be
%    consulted, for example Peerlings (lecture notes 4A700).
%
%    Local number of a numbering of 2-D elements:
%
%    Linear triangular element:
%          3   
%         / \
%        /   \ 
%       /     \
%      1-------2
%    then CONN = [1 2 3].
%
%    Quadratic triangular element:
%          3   
%         / \
%        6   5 
%       /     \
%      1---4---2
%    then CONN = [1 2 3 4 5 6].    
%
%    Linear quadrilateral element:
%      4-------3
%      |       |
%      |       |
%      |       |
%      1-------2
%    then CONN = [1 2 3 4].
%
%    Quadratic quadrilaterals element, Serendipity:
%      4---7---3
%      |       |
%      8       6
%      |       |
%      1---5---2
%    then CONN = [1 2 3 4 5 6 7 8].
%
%    Quadratic quadrilaterals element, Lagrange:
%      4---7---3
%      |       |
%      8   9   6
%      |       |
%      1---5---2
%    then CONN = [1 2 3 4 5 6 7 8 9].
%  
%    [H,HM,HN,HE] = FEMPLOT( ... , 'PropertyName' , property) different
%    options can be selected:
%  
%    NODES
%      {'on'} | 'off'
%    show marker at the nodal positions.
%  
%    NODESNUMBERS
%      {'on'} | 'off'%    display the node numbers
%
%    ELEMENTNUMBERS
%      {'on'} | 'off'%    display the element numbers
%
%    COLOR
%      select a different color for the element edges. By default the
%      color is blue.
%
%    Secondly, the different handles of the mesh-plot can be returned. H
%    is the handle of the patches. HM is the handles of the markers
%    at the nodes that do not coincide with element edges (optional:
%    depends on the type of elements). HN and HE are the handles of the
%    node and  element numbers respectively (optional: depends on the
%    option selected).
% 
% SEE ALSO: TENSOR/FEMGRID 

% Set defaults
nds = 0;   % markers at node positions
ndn = 0;   % node numbers
eln = 0;   % element numbers
fc  = 'b'; % facecolor
h   = [];  % handle of the patch
hm  = [];  % handle of the mid-point nodes
hn  = [];  % handle of node numbers
he  = [];  % handle of element numbers

nar=size(varargin,2); i=0;
while i<nar, i=i+1;
    if ~ischar(varargin{i}), error('Unknown option'); 
    elseif strcmpi(varargin{i},'Nodes'); i=i+1;
        if ~ischar(varargin{i}), error('Unknown NODES option');
        elseif strcmpi(varargin{i},'on'), nds=1;
        elseif strcmpi(varargin{i},'off'), nds=0;
        else, error('Unknown NODES option'); 
        end
    elseif strcmpi(varargin{i},'NodeNumbers'); i=i+1;
        if ~ischar(varargin{i}), error('Unknown NODENUMBERS option');
        elseif strcmpi(varargin{i},'on'), ndn=1;
        elseif strcmpi(varargin{i},'off'), ndn=0;
        else, error('Unknown NODES option');
        end
    elseif strcmpi(varargin{i},'ElementNumbers'); i=i+1;
        if ~ischar(varargin{i}), error('Unknown ELEMENTNUMBERS option');
        elseif strcmpi(varargin{i},'on'), eln=1;
        elseif strcmpi(varargin{i},'off'), eln=0;
        else, error('Unknown NODES option');
        end
    elseif strcmpi(varargin{i},'Color'); i=i+1; fc=varargin{i};
    else, error('Unkown option');
    end
end

% read input
dim = size(x.basis,1);
t0  = t; % backup connectivity
p   = reshape(x.components, size(x.components,1)/size(x.basis,1), size(x.basis,1) ); % nodal positions in doubles

% reshape t: elements are split to faces

if dim == 2
	if size(t,2)==3 % 2D: linear triangles
		tm=[]; 
		t =t;
	elseif size(t,2)==6 % 2D: bi-quadratic triangles
		tm=[];
		t =t(:,[1 4 2 5 3 6]);
	elseif size(t,2)==4 % 2D: linear quadrilaterals
		tm=[];
		t =t;
	elseif size(t,2)==8 % 2D: Serendipity: bi-quadratic quadrilaterals
		tm=[];
		t =t(:,[1 5 2 6 3 7 4 8 1]);
	elseif size(t,2)==9 % 2D: Lagrange: quadratic quadrilaterals
		tm=t(:,9);
		t =t(:,[1 5 2 6 3 7 4 8 1]);
	else, error('Undefined element');
	end
    
elseif dim == 3
	if size(t,2)==4 % 3D: linear triangular elements
		tm=[];
		t =[t(:,[1 2 3]);
            t(:,[1 4 2]);
            t(:,[4 3 2]);
            t(:,[1 3 4])];
	elseif size(t,2)==10 % 3D: bi-quadratic triangular elements
		tm=[];
		t =[t(:,[1 5 2 9 4 8]);
            t(:,[2 6 3 10 4 9]);
            t(:,[1 7 3 10 4 8]);
            t(:,[1 5 2 6 3 7])];
	elseif size(t,2)==8 % 3D: linear quadrilaterals
		tm=[];
		t =[t(:,[1 2 3 4]);
            t(:,[1 2 6 5]);
	 		t(:,[2 6 7 3]);
 			t(:,[3 7 8 4]);
			t(:,[4 8 5 1]);
			t(:,[5 6 7 8]) ];
	elseif size(t,2)==20 % 3D: Seredipity: bi-quadratic quadrilaterals
		tm=[];
		t =[t(:,[ 1  9  2 10  3 11  4 12]);
			t(:,[ 5 13  6 14  7 15  8 16]);
			t(:,[ 1  9  2 18  6 13  5 17]);
			t(:,[ 2 10  3 19  7 14  6 18]);
			t(:,[ 3 11  4 20  8 15  7 19]);
			t(:,[ 4 20  8 16  5 17  1 12]) ];
	elseif size(t,2)==27 % 3D: Lagrange: quadratic quadrilaterals
		tm= t(:,[21 22 23 24 25 26 27   ]);
		t =[t(:,[ 1  9  2 10  3 11  4 12]);
			t(:,[ 5 13  6 14  7 15  8 16]);
			t(:,[ 1  9  2 18  6 13  5 17]);
			t(:,[ 2 10  3 19  7 14  6 18]);
			t(:,[ 3 11  4 20  8 15  7 19]);
			t(:,[ 4 20  8 16  5 17  1 12]) ];
	else, error('Undefined element');
	end
	
else, error('Incorrect dimension')
end

% sweep t and tm
tmp=sort(t,2); [tmp,in]=unique(tmp,'rows'); t=t(in,:);
tm=tm(:); tm=unique(tm);


% plot element
if ~nds
    h=patch('Faces',t,'Vertices',x,'EdgeColor',fc,'FaceColor','none');
else
	hold on;
    h=patch('Faces',t,'Vertices',x,'EdgeColor',fc,'FaceColor','none','Marker','o','MarkerEdgeColor',fc,'MarkerFaceColor',fc);
    if     dim==2, hm=plot (p(tm,1),p(tm,2)        ,'Marker','o','MarkerEdgeColor',fc,'MarkerFaceColor',fc,'linestyle','none');
	elseif dim==3, hm=plot3(p(tm,1),p(tm,2),p(tm,3),'Marker','o','MarkerEdgeColor',fc,'MarkerFaceColor',fc,'linestyle','none');
	end
end

% for a 3D plot, reset the view
if dim == 3
    view(45,45);
end

% Plot the node numbers

if ndn
	hold on;

    % get a list with node numbers
    nn=t0(:); nn=unique(nn);
	tx=cell(length(nn),1);
    for i=1:length(nn), tx{i}=sprintf('%d',nn(i)); end
	
	if dim==2,     hn=text(p(nn,1),p(nn,2)       ,tx);
	elseif dim==3, hn=text(p(nn,1),p(nn,2),p(nn,3),tx);
	end

	set(hn,'Color',fc);

end

if eln
	hold on;

	pav=zeros(size(t0,1),dim);
	tx=cell(size(t0,1),1);
	for i=1:size(t0,1), pav(i,:)=mean(p(t0(i,:),:),1); tx{i}=sprintf('%d',i); end

	if dim==2,     he=text(pav(:,1),pav(:,2)         ,tx);
	elseif dim==3, he=text(pav(:,1),pav(:,2),pav(:,3),tx);
	end

    set(he,'Color',fc,'EdgeColor',fc);

end

varargout{1}=h;  % face handle
varargout{2}=hm; % handle of midpoint nodes
varargout{3}=hn; % handle of node numbers
varargout{4}=he; % handle of element numbers

if nargout>4, error('Too many output arguments'); end

end
