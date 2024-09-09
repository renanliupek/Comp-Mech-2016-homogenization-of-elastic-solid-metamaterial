%% Renan Liupekevicius Carnielli TU/e
% start in 28-04-2022
% computational homogenization of LRAM
clear; close all;


% TOPOLOGY: no Fluid and solid a must at the corners
% TRANSPARENT ELEMENTS: DOES NOT HANDLE transparent massless solid elements




%% IMPORT MESH
disp('codeblock: SELECT DESIGN')
% select unit cell design
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design = 'LRAM_coarse';
disp(design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------



% read mesh text file comsol 5.4
%--------------------------------------------------------------------------
disp('codeblock: READ & CONVERT LINEAR TO QUADRATIC ELEMENT MESH')
% read text file
l   = read_mphtxt_54(design);

% convert to quadratic element
mesh=mesh2d_lin2qua_uc(l);
%--------------------------------------------------------------------------


% read mesh text file comsol 5.6
%--------------------------------------------------------------------------
% disp('codeblock: READ MESH COMSOL 5.6')
% mesh   = read_mphtxt_56(design);
%--------------------------------------------------------------------------

% copy to local variables (please double click on 'mesh' struct in 
% workspace to see the meaning of each cell)
%--------------------------------------------------------------------------
  x    = mesh.x{2,2};    % tensor form
  mx   = mesh.x{2,1};    % matrix form
  conn = mesh.elem{2,1}; % conn matrix of element
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% mesh parameters
  m  = size(  conn, 1); % number of elements in the mesh
  n  = size(  x, 1);    % number of nodes in the mesh 
  fprintf('NUMBER OF MESH ELEMENTS: %d\n', m)
  fprintf('NUMBER OF MESH NODES: %d\n', n)
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit cell size
  a = 2*max(  mx(:,1));
ax=a;
ay=a;

if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not x centralized');end
if(max(mx(:,1))+min(mx(:,1))>1e-8); warning('RVE is not y centralized');end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

%% INITIAL DEFINITIONS MATERIAL PARAMETERS

%--------------------------------------------------------------------------
dim   = 2;  % problem dimension
thick = 1;  % thickness (out-of-plane direction), in [m]
%--------------------------------------------------------------------------






% IF Z. LIU (SCIENCE 2000) DESIGN
if strncmp(design,'LRAM',4) % LRAM's mesh

    n_solid_phases = 3;
    % define solid phase parameters
    %----------------------------------------------------------------------
    %matrix (part 1)
    matProp{1}.E       =   3.6e9; % Young's modulus
    % matProp{1}.E       =   3.6e6; warning(' test soft matrix')
    matProp{1}.nu      =   0.368; % poisson ratio
    matProp{1}.rho     =   1180; % mass density
    matProp{1}.G       =   matProp{1}.E/2/(1 +   matProp{1}.nu);
    matProp{1}.kappa   =   matProp{1}.E/3/(1 - 2*matProp{1}.nu);
    matProp{1}.cL      =  sqrt(matProp{1}.E * (1-matProp{1}.nu) ...
                                /(1+matProp{1}.nu)/(1 - 2*matProp{1}.nu)...
                                /matProp{1}.rho);
    
    
    % coating (part 2)
    matProp{2}.E       =   0.1175e6; % original
    matProp{2}.nu      =   0.469; %MAKE IT nu=0.3 TO ENSURE NON-COMPRESSIB.
    matProp{2}.rho     =   1300;
    matProp{2}.G       =   matProp{2}.E/2/(1 +   matProp{2}.nu);
    matProp{2}.kappa   =   matProp{2}.E/3/(1 - 2*matProp{2}.nu);
    matProp{2}.cL      =  sqrt(matProp{2}.E * (1-matProp{2}.nu) ...
                                /(1+matProp{2}.nu)/(1 - 2*matProp{2}.nu)...
                                /matProp{2}.rho);
    
    % core (part 3)
    matProp{3}.E       =   40.82e9;
    matProp{3}.nu      =   0.37;
    matProp{3}.rho     =   11600;
    matProp{3}.G       =   matProp{3}.E/2/(1 +   matProp{3}.nu);
    matProp{3}.kappa   =   matProp{3}.E/3/(1 - 2*matProp{3}.nu);
    matProp{3}.cL      =  sqrt(matProp{3}.E * (1-matProp{3}.nu) ...
                                /(1+matProp{3}.nu)/(1 - 2*matProp{3}.nu)...
                                /matProp{3}.rho);
    %----------------------------------------------------------------------


% IF METAFOAM DESIGN
    elseif strncmp(design,'design6',7)
    
    n_solid_phases = 2;
    % define solid phase parameters
    %----------------------------------------------------------------------
    %matrix (part 1)
    matProp{1}.E       =   1e6; % Young's modulus
    matProp{1}.nu      =   0.4; % poisson ratio
    matProp{1}.rho     =   1000; % mass density
    matProp{1}.G       =   matProp{1}.E/2/(1 +   matProp{1}.nu);
    matProp{1}.kappa   =   matProp{1}.E/3/(1 - 2*matProp{1}.nu);
    
    
    % coating (part 2)
    matProp{2}.E       =   180e9; % original
    % matProp{2}.E       =   3e6; % test light mass
    matProp{2}.nu      =   0.3; %
    matProp{2}.rho     =   1e6;
    % matProp{2}.rho     =   3000; % test light mass
    matProp{2}.G       =   matProp{2}.E/2/(1 +   matProp{2}.nu);
    matProp{2}.kappa   =   matProp{2}.E/3/(1 - 2*matProp{2}.nu);
    %----------------------------------------------------------------------
   
    % string stiffness
    %----------------------------------------------------------------------
    E = matProp{1}.E;
    %----------------------------------------------------------------------


% ELSE: CONVENTIONAL FOAM (1 SOLID PHASE)
else
    %--------------------------------------------------------------------------
    %solid material properties PU from Mirka ch3 (IN TERMS OF E AND nu)
      rhos     = 1000;   % in [kg/m3]
      E        = 1e6;    % [Pa]
      nu       = 0.4;    % []    
      G        =   E/2/(1+    nu);
      kappa    =   E/3/(1-2*  nu);
    %--------------------------------------------------------------------------
    
end





% define tensor basis
%--------------------------------------------------------------------------
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(  b{1},   b{2});
e1 = ee(1);
e2 = ee(2);
%--------------------------------------------------------------------------


% calculate fundamental tensors
%--------------------------------------------------------------------------
  I   = identity(2, ee);
  I4  = identity(4, ee);
  I4S = 1/2 * (  I4 + rtranspose( I4));
%--------------------------------------------------------------------------


%%  DESIGNs: LIST OF TAGS 

% MANUALLY import tags from comsol: check acoustic-structure boundary
% list under 'boundary selection'. Same procedure for selecting tag-list of
% the solid phase(s).


switch sscanf(design, '%c')
    
    case 'LRAM' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'LRAM_coarse' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'LRAM_coarser' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------


    case 'LRAM_coarse_shifted' 
    % same as LRAM_coarse shifted (a/2,a/2)
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [3 4 7 8];  
    % solid elements tag
    list_solid2_tag = [2 5 9 10];  
    % solid elements tag
    list_solid3_tag = [1 6 11 12];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add here new case if you designed a new geometry/mesh
% case 'new_mesh_here' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% SORTING SOLID 

% TOTAL SOLID ELEMENTS DISREGARDING FROM WHICH MATERIAL PHASE
%--------------------------------------------------------------------------    
    % solid elements
      solid_elems =[];
      for i=list_solid_tag
      solid_elems = [solid_elems find( mesh.elem{2,2}==i).'];
      end
    % position vector of the solid nodes
      conns       =   conn(  solid_elems,:);
    % solid nodes
      nodes_s     = unique(reshape(  conns,1,[]));
    % number of solid elements in the mesh
      ms          = size(  conns,1);
    % number of solid nodes
      ns          = length(  nodes_s);
%--------------------------------------------------------------------------


% if LRAM's mesh
if strncmp(design,'LRAM',4)

    % SOLID 1
    %--------------------------------------------------------------------------    
        % solid elements
          solid1_elems =[];
          for i=list_solid1_tag
          solid1_elems = [solid1_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{1}       =   conn(  solid1_elems,:);
        % solid nodes
          nodes_ss{1}     = unique(reshape(  connss{1},1,[]));
        % number of solid elements in the mesh
          mss{1}          = size  (  connss{1}, 1);
        % number of solid nodes
          nss{1}          = length(  nodes_ss{1});
    %--------------------------------------------------------------------------
    
    % SOLID 2
    %--------------------------------------------------------------------------    
        % solid elements
          solid2_elems =[];
          for i=list_solid2_tag
          solid2_elems = [solid2_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{2}       =   conn(  solid2_elems,:);
        % solid nodes
          nodes_ss{2}     = unique(reshape(  connss{2},1,[]));
        % number of solid elements in the mesh
          mss{2}          = size  (  connss{2}, 1);
        % number of solid nodes
          nss{2}          = length(  nodes_ss{2});
    %--------------------------------------------------------------------------
    
    % SOLID 3
    %--------------------------------------------------------------------------    
        % solid elements
          solid3_elems =[];
          for i=list_solid3_tag
          solid3_elems = [solid3_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{3}       =   conn(  solid3_elems,:);
        % solid nodes
          nodes_ss{3}     = unique(reshape(  connss{3},1,[]));
        % number of solid elements in the mesh
          mss{3}          = size  (  connss{3},1);
        % number of solid nodes
          nss{3     }     = length(  nodes_ss{3});
    %--------------------------------------------------------------------------

% if metafoam mesh    
elseif strcmp(design,'design6') 

    
    % SOLID 1
    %--------------------------------------------------------------------------    
        % solid elements
          solid1_elems =[];
          for i=list_solid1_tag
          solid1_elems = [solid1_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{1}       =   conn(  solid1_elems,:);
        % solid nodes
          nodes_ss{1}     = unique(reshape(  connss{1},1,[]));
        % number of solid elements in the mesh
          mss{1}          = size  (  connss{1}, 1);
        % number of solid nodes
          nss{1}          = length(  nodes_ss{1});
    %--------------------------------------------------------------------------
    
    % SOLID 2
    %--------------------------------------------------------------------------    
        % solid elements
          solid2_elems =[];
          for i=list_solid2_tag
          solid2_elems = [solid2_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{2}       =   conn(  solid2_elems,:);
        % solid nodes
          nodes_ss{2}     = unique(reshape(  connss{2},1,[]));
        % number of solid elements in the mesh
          mss{2}          = size  (  connss{2}, 1);
        % number of solid nodes
          nss{2}          = length(  nodes_ss{2});
    %--------------------------------------------------------------------------


end


% FLUID
%--------------------------------------------------------------------------
    % fluid elements
      fluid_elems = setdiff(1:  m,  solid_elems);
    % conn fluid
      connf       =   conn(  fluid_elems,:);
    % fluid nodes
      nodes_f     = unique(reshape(  connf,1,[]));
    % number of fluid elements in the mesh
      mf          = size(  connf,1);
    % number of fluid nodes
      nf          = length(  nodes_f);
%--------------------------------------------------------------------------



%% PLOT MESH

% if LRAM's mesh
if strncmp(design,'LRAM',4)

    %--------------------------------------------------------------------------
    % plot each solid phase
    figure(2)
    clf
    hold on;
    femplot(  x,  connss{1},'Color' , 'black', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    femplot(  x,  connss{2},'Color' , 'red', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    femplot(  x,  connss{3},'Color' , 'blue', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % scatter(mx(1:2,1),mx(1:2,2))
    hold off;
    % exportgraphics(gca,'plot.png','BackgroundColor','none')
    %--------------------------------------------------------------------------




% if metafoam mesh    
elseif strncmp(design,'design6',7)

    %--------------------------------------------------------------------------
    % plot each solid phase and fluid phase
    figure(2)
    clf
    hold on;
    femplot(  x,  connf,'Color' , 'blue', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    femplot(  x,  connss{1},'Color' , 'black', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    femplot(  x,  connss{2},'Color' , 'red', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % scatter(mx(1:2,1),mx(1:2,2))
    hold off;
    % exportgraphics(gca,'plot.png','BackgroundColor','none')
    %--------------------------------------------------------------------------



% else conventional foam mesh
else

    %--------------------------------------------------------------------------
    %plot each phase
    figure(2)
    clf
    hold on;
    femplot(  x,  conns,'Color' , 'black', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    femplot(  x,  connf,'Color' , 'blue', 'Nodes', 'off', ...
    'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    hold off;
    % exportgraphics(gca,'plot.png','BackgroundColor','none')
    %--------------------------------------------------------------------------

end






%% NODE SORTING OF BOUNDARIES 
disp('codeblock: NODE SORTING OF BOUNDARIES ')

% NOTE: this codeblock assumes the top and bottom boundaries are solid
% boundaries.

%--------------------------------------------------------------------------
% mesh precision. Check in COMSOL for correct use of it. Here I set 1e-8 by
% try and error which is not ideal.
precision = 1e-8; 
%--------------------------------------------------------------------------


% sort nodes of the solid and fluid boundaries
%--------------------------------------------------------------------------
% nodes on the left boundary of the mesh: line x = -a/2
[  l,~]        = find(abs(  mx(:,1) +   a/2 )<precision ); % 
  left_s       = intersect(  l,  nodes_s)';  % solid boundary

% nodes on the right boundary of the mesh: line x = +a/2                  
[  r,~]        = find(abs(  mx(:,1) -   a/2 )<precision ); % 
  right_s      = intersect(  r,  nodes_s)';  % fluid boundary

% nodes on the bottom boundary of the mesh: line y = -a/2
[  bottom_s,~] = find(abs(  mx(:,2) +   a/2 )<precision ); % 

% nodes on the bottom boundary of the mesh: line y = +a/2
[  top_s,~]    = find(abs(  mx(:,2) -   a/2 )<precision ); % 


% corner nodes
  corners_s = [   intersect(  l,  bottom_s) ...
                  intersect(  r,  bottom_s) ...
                  intersect(  r,  top_s   ) ...
                  intersect(  l,  top_s   )];
% prescribed corners
  pcorners = [  corners_s(1)   corners_s(2)   corners_s(4)];


% exclude corners from solid node lists
  left_s   = setdiff(  left_s  ,   corners_s);
  right_s  = setdiff(  right_s ,   corners_s);
  top_s    = setdiff(  top_s   ,   corners_s)';
  bottom_s = setdiff(  bottom_s,   corners_s)';

clear   l; clear   r;
%--------------------------------------------------------------------------


%% INDEX TABLE
% __________________________________________________________________________
% FLUID (F)
% SOLID (S)
%
%        CURRENT_INDEX      
%               |  
%               V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________


%%  4 INTEGRATION POINTS WITHIN MASTER ELEMENT
%--------------------------------------------------------------------------
% Notation from TensorLab
% quadratic quadrilateral element:
%         4----7----3
%         |         |
%         8         6
%         |         |
%         1----5----2
%      then CONN = [1 2 3 4 5 6 7 8].
%
% The coordinates of integration points x (in the coordinate system of the
% master element). Exact integration up to polynomial of order 5.
%          --------- 
%         | x  x  x |
%         | x  x  x |
%         | x  x  x |
%          ---------
% and the weight factors
%--------------------------------------------------------------------------
xi = [ -sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 + sqrt(3/5) *e2
       -sqrt(3/5) *e1 + sqrt(3/5) *e2
                0 *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 +         0 *e2       
                0 *e1 + sqrt(3/5) *e2                
       -sqrt(3/5) *e1 +         0 *e2
                0 *e1 +         0 *e2   ];
w  = [ 25/81
       25/81
       25/81
       25/81
       40/81
       40/81
       40/81
       40/81
       64/81     ];
%--------------------------------------------------------------------------


%% ASSEMBLE STIFFNESS AND MASS MATRICES (SOLID ELEMENTS)
disp('codeblock: ASSEMBLING SOLID MASS AND STIFFNESS MATRICES')
%--------------------------------------------------------------------------
% Knn and Mnn are only temporary matrix and they are not the matrices K and
% M from discretized form in the notes. Note their size is the total number
% of nodes which include the fluid nodes. It's only relevant property is
% that it preserves the node indexing.
%--------------------------------------------------------------------------


% declare stiffness matrix tensor with initially zero-2nd-order tensors
%--------------------------------------------------------------------------
Knn = zeros(2, ee,   n,   n); 
Mnn = zeros(  n,   n);
%--------------------------------------------------------------------------



% SPRING CONNECTIVITY 
%--------------------------------------------------------------------------
% disp('ASSEMBLED SPRING ELEMENTS'); 
%     for e = 1:size(  connspring,1) % loop over all string elements
%     
%     % extract node position from current string element   
%     iie  =   connspring(e, :); % nodes of the curren  t element e
%     xe   =   x(iie);      % coordinates of these nodes
%     
%     % unit vector in i-j direction in the undeformed configuration
%     l0ij = (xe(2)-xe(1))/norm(xe(2)-xe(1));
% 
%     % compute elementary matrices (maybe minus)                   
%     Ke =  ks*l0ij*l0ij* [    1    -1 
%                             -1     1 ];
% 
%     % assembly
%     Knn(iie, iie) = Knn(iie, iie) + Ke;
% 
%    end % end of element loop
%--------------------------------------------------------------------------




% if multiple solid phases
if strncmp(design,'LRAM',4)
% SOLID ELEMENTS ( MULTIPLE PHASE MODEL)
%--------------------------------------------------------------------------

%initialize porosity to 1
phi = 1;

%loop
for i =1:n_solid_phases

 % define material parameters from phase i   
 tC    =  matProp{i}.kappa*I*I + 2* matProp{i}.G*(  I4S - 1/3*  I*  I);
 rhos  =  matProp{i}.rho;

 % initialize area of each solid phase to zero
 V{i}=0;

    for e = 1:mss{i} % loop over all solid elements of phase i
    
    % % display computation percentage
    % if mod(floor(e/mss{i}*100),10)==0
    % fprintf('solid %d: assembled %.2f ', i, e/mss{i}*100); disp('%');
    % end
    
    % extract nodes
    iie =   connss{i}(e, :); % nodes of the current element e
    xe  =   x(iie);            % coordinates of these nodes
        
    % compute element matrices in integration point loop
    Ke = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee.
    Me = zeros(8,8);         % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ke = Ke + w(k) * dot(gradNe,   tC, gradNe') *   thick * det(J);
        Me = Me + w(k) *  Ne *  rhos*   Ne'  *   thick * det(J);     
   
    end % end of interation point loop
   
   % assembly
    Knn(iie, iie) = Knn(iie, iie) + Ke;
    Mnn(iie, iie) = Mnn(iie, iie) + Me;


   % get element area
   V{i}=V{i}+get_element_area(xe,ee);

   end % end of element loop

% compute porosity   
phi = phi - V{i}/a^2;

% assignment actual mass and stiffness matrices
K=Knn(  nodes_s,  nodes_s);
M=Mnn(  nodes_s,  nodes_s);
end

%--------------------------------------------------------------------------


else


% SOLID ELEMENTS ( ONE SOLID PHASE MODEL)
%--------------------------------------------------------------------------
  % initialize area of solid phase to zero
  V=0;

  % define stiffness tensor
  tC   =   kappa *  I*  I + 2*  G * (  I4S - 1/3*  I*  I);

  for e = 1:size(  conns,1) % loop over all solid elements
    
    % display computation percentage
    if mod(floor(e/ms*100),10)==0
    fprintf('solid: assembled %.2f ', e/ms*100);disp('%');
    end


    iie =   conns(e, :); % nodes of the current element e
    xe  =   x(iie);     % coordinates of these nodes
        
    % compute element matrices in integration point loop
    Ke = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee.
    Me = zeros(8,8);         % zero matrix of 8x8 scalars

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Ke = Ke + w(k) * dot(gradNe,   tC, gradNe') *   thick * det(J);
        Me = Me + w(k) *  Ne *  rhos*   Ne'  *   thick * det(J);     
        


    end % end of interation point loop
   

   % assembly
    Knn(iie, iie) = Knn(iie, iie) + Ke;
    Mnn(iie, iie) = Mnn(iie, iie) + Me;

   % get element area
   V=V+get_element_area(xe,ee);    
   

   end % end of element loop

% assignment actual mass and stiffness matrices
K=Knn(  nodes_s,  nodes_s);
M=Mnn(  nodes_s,  nodes_s);

% compute porosity
phi=(a^2-V)/a^2;
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% force symmetrization ( "correct" numerical errors)
K = (K+K')/2;
M = (M+M')/2;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
clear Me; clear Ke;
clear Knn; clear Mnn;
%--------------------------------------------------------------------------







%% PRESCRIBED AND FREE NODES SPLIT
%--------------------------------------------------------------------------
% INDEX OF DISPLACEMENT DOF IN NODES_S (ELIMINATED DEPENDENT DOFS)
% same index change procedure must be executed but only in the solid phase.

%|__NOTES_NOTATION_|____ DESCRIPTION_______________________________________
%|                 |
%|        p        | prescribed nodes in the solid phase
%|_________________|_______________________________________________________
%|                 |
%|        f        | free nodes in the solid phase                                                      
%|                 |
%|_________________|_______________________________________________________
%|                 |
%|        p'       | prescribed nodes in the fluid phase
%|_________________|_______________________________________________________
%|                 |
%|        f'       | free nodes in the fluid phase
%|_________________|_______________________________________________________
%--------------------------------------------------------------------------


% %--------------------------------------------------------------------------
% % solid
% p  = get_ind(dofs_re,dofs_pc); % = [1 2 3] by construction
% f  = setdiff(1:ns_re,p);
% 
% % fluid
% pp = [doff_l doff_r];
% fp = setdiff(1:  nf,pp);
% %--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% solid
p  = [corners_s left_s right_s bottom_s top_s]; % = [1 2 3] by construction
f  = setdiff(1:ns,p);

%--------------------------------------------------------------------------


%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
% SOLID (S)
%
%                                                            CURRENT_INDEX      
%                                                                   |  
%                                                                   V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   | nodes_f |   nodes_f    |    nodes_f   |  nodes_f
%____________|__(FS) _|_(S)_(F)_|__(S)_____(F)_|___(S)___(F)__|__(S)_(F)__
%            |    1   |   1 1   |  dofs_c   1  |  dofs_pc  1  |    p  pp
%            |    2   |   2 2   |  dofs_un  2  |  dofs_un  2  |    f  ff
%   ORDER    |    .   |   . .   |  dofs_in  .  |  dofs_in  .  |
%            |    .   |   . .   |  dofs_de  .  |           .  |
%            |    .   |   . .   |           nf |           nf | 
%            |    n   |  ns nf  |              |              |
%____________|___ ____|_________|______________|______________|____________

%% PARTITIONING

%--------------------------------------------------------------------------
% solid
 
    
    % solid
    % mass tensor matrix
    % M =  Mpp Mpf   
    %      Mfp Mff 
    M_p_p   = M(p , p ); M_p_f   = M(p , f );
    M_f_p   = M(f , p ); M_f_f   = M(f , f );
    % stiffness tensor matrix
    % K =  Kpp Kpf   
    %      Kfp Kff 
    K_p_p   = K(p , p ); K_p_f   = K(p , f );
    K_f_p   = K(f , p ); K_f_f   = K(f , f );


%--------------------------------------------------------------------------



%% TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX
disp('codeblock: TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX')

%--------------------------------------------------------------------------

% mass matrix M
sM_p_p = tm2sm_v( M_p_p*I, ee , 2 ); sM_p_f = tm2sm_v( M_p_f*I, ee , 2 );
sM_f_p = tm2sm_v( M_f_p*I, ee , 2 ); sM_f_f = tm2sm_v( M_f_f*I, ee , 2 );

% stiffness matrix K
sK_p_p = tm2sm_v( K_p_p, ee , 2 ); sK_p_f = tm2sm_v( K_p_f, ee , 2 );
sK_f_p = tm2sm_v( K_f_p, ee , 2 ); sK_f_f = tm2sm_v( K_f_f, ee , 2 );

%--------------------------------------------------------------------------


%% ASSEMBLE SYMMETRIC SYSTEM OF EQUATIONS 
disp('codeblock: SOLVE')

% recap on -C^T relations
%--------------------------------------------------------------------------
%     C_pp_p = - C_p_pp.' ; C_pp_f = - C_f_pp.' ;
%     C_fp_p = - C_p_fp.' ; C_fp_f = - C_f_fp.' ;
%--------------------------------------------------------------------------


% lamda submatrices
%--------------------------------------------------------------------------
lam_P_P = [  sM_p_p ];

lam_P_F = [  sM_p_f ];

lam_F_P = [  sM_f_p ];

lam_F_F = [  sM_f_f ];
%--------------------------------------------------------------------------


% mu submatrices
%--------------------------------------------------------------------------
mu_P_P  = [  sK_p_p ];  

mu_P_F  = [  sK_p_f ];  

mu_F_P  = [  sK_f_p ];  

mu_F_F  = [  sK_f_f ];  
%--------------------------------------------------------------------------



% warning('symmetrization of K- and M- blocks forced to get rid of numerical error')
% % force symmetrization if only-solid problem or only-fluid problem
% %--------------------------------------------------------------------------
%     disp('SYMMETRIZATION FORCED')
%     mu_F_F   = sparse((mu_F_F + mu_F_F.' )/2);
%     lam_F_F  = sparse((lam_F_F+ lam_F_F.')/2);
%     mu_P_P   = sparse((mu_P_P + mu_P_P.' )/2);
%     lam_P_P  = sparse((lam_P_P+ lam_P_P.')/2);
%     mu_P_F   = mu_F_P.';
%     lam_P_F  = lam_F_P.';
% %--------------------------------------------------------------------------




%% EIGEN: EIGEN STUDY, PRINT EIGENFREQUENCIES AND EIGENMODES

% norm(mu_F_F, 'fro')
% norm(mu_F_F - sparse((mu_F_F + mu_F_F.')/2), 'fro')
% norm(lam_F_F, 'fro')
% norm(lam_F_F- sparse((lam_F_F+ lam_F_F.')/2), 'fro')


%--------------------------------------------------------------------------
% define desired number of modes to compute/display
n_modes = 6;
disp('ITERATIVE SOLVER')
    % Subset of eigenvalues -> CHEAPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % eigs is recommended by matlab when sparse nonsymmetric matrices.
    % % compute eigen modes V, and eigen values D with smallest absolute
    % value right problem.
    [phi_F_Q, Dr] = eigs(mu_F_F ,lam_F_F ,n_modes,...
                   'smallestabs',...
                   'Display',0,...
                   'FailureTreatment','replacenan');
    % left problem (psi appears always transposed along the expressions)
%     [psi_F_Q, Dl] = eigs(mu_F_F',lam_F_F',n_modes,...
%                    'smallestabs',...
%                    'Display',0,...
%                    'FailureTreatment','replacenan');
    psi_F_Q = phi_F_Q;  
  

    % CHECK whether the frequencies in Dr and Dl are enoughly equal. For
    % the current example they are equal up to the 5th digit after comma.
    % log10(abs(Dr-Dl))

    % print frequencies f
    freqs = diag(sqrt(Dr)/(2*pi));
    fprintf('Eigen frequencies:\n');
    fprintf('%f Hz\n', freqs);
%--------------------------------------------------------------------------


%direct solver
%--------------------------------------------------------------------------
% disp('DIRECT SOLVER')
% [phi_F_Q, Dr] = eig(full(mu_F_F) ,full(lam_F_F));
% % print frequencies f
% freqs = sort(diag(sqrt(Dr)/(2*pi)));
% fprintf('Eigen frequencies:\n');
% fprintf('%f Hz\n', freqs(1:n_modes));
%--------------------------------------------------------------------------
    

% EIGEN:  PHASE CORRECTION
%--------------------------------------------------------------------------
% IMPORTANT NOTE: PHASE CORRECTION BEFORE NORMALIZATION 
% Compute round(psi'*lamFF*phi). It should be positive diag matrix - note 
% that this product should be approx positive diag because the the
% frequencies must be positive-. Since the right and left eigen modes were
% computed with different calls from the function eigs, they can
% potentially converge to a shape mode with phase difference of pi. In 
% this case, psi'*lamFF*phi ( or round(psi'*lamFF*phi) for getting rid 
% of almost zero off diag elements) may have some of the diagonal elements
% as -negative instead of +positive. To correct this effect it is necessary
% to multiply either phi or psi by the 'sign' matrix 'I'= psi'*lamFF*phi.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% get diagonal
diag_psi_lam_phi   = diag(psi_F_Q'*lam_F_F*phi_F_Q);
% make it a vector with only the sign of each element
diag_psi_lam_phi   = diag_psi_lam_phi./abs(diag_psi_lam_phi);
% make 'identity' matrix with sign correction
I_phase_correction = diag(diag_psi_lam_phi);

% correct the sign of left eigen modes with respect to the right modes
psi_F_Q = psi_F_Q * I_phase_correction;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% save version of phi for plotting
phi_plot = phi_F_Q;
%--------------------------------------------------------------------------

% EIGEN:  MODE NORMALIZATION WITH RESPESCT TO THE LAMDA MATRIX
%--------------------------------------------------------------------------
% ready to mass normalization, now, the product psi'*lamFF*phi has only
% postive diagonal compared to previous case
% get norm of each mode


vec_norms = sqrt(diag(psi_F_Q'*lam_F_F*phi_F_Q));

% figure(8);subplot(2,1,1); loglog(abs(vec_norms)); hold on;

phi_F_Q=phi_F_Q./vec_norms'; %divide each column of phi by element of vec_norms 
psi_F_Q=psi_F_Q./vec_norms'; %divide each column of psi by element of vec_norms 

% subplot(2,1,2); loglog(abs(sqrt(diag(psi_F_Q'*lam_F_F*phi_F_Q))));

% CONSISTENCY CHECK: 
% Compute again psi'*lamFF*phi. It should be approx. the identity matrix.
% when eigs is used the phase problem doen't exist. It can be seen by
% computing psi(:,1:n_modes)'*lamFF*phi(:,1:n_modes). The diagonal is made
% out of positive integers.
%--------------------------------------------------------------------------





% EIGEN: ASSIGN DISP. SHAPE MODES TO TENSOR VECTOR & PLOT
%--------------------------------------------------------------------------
% displacement 'ured' is followed is 'u' displacement reduced 'red' thanks
% to the elimitation of the dependend points due to periodic bc
% 'ured' stands for 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% prescribed (reduced) displacements are zero
    % declare & assign
su_p = zeros( 2*length(p), size(phi_F_Q,2) ); 
        
% free (reduced) displacements are the eigenvectors    
% u_f = phi_plot(1:2:2*length(f)-1,:)*e1+ phi_plot(2:2:2*length(f),:)*e2;
su_f = phi_plot;


su_pf_Q = [su_p
           su_f];

vu_pf_Q = su_pf_Q(1:2:end,:)*e1+ su_pf_Q(2:2:end,:)*e2;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% prepare for plot: 
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot = zeros(1,ee,length(p), size(phi_F_Q,2));
    vu_plot([p f],:) = vu_pf_Q;
%--------------------------------------------------------------------------







% PLOT EIGEN DISP MODE SHAPES (RESCALED DISP.)
%--------------------------------------------------------------------------
disp('codeblock: PLOT')



sfreqs=freqs; %sorted frequencies vector
mode=6;




% if LRAM's geometry
if strncmp(design,'LRAM',4)


% plot the solid displacement in case only 3 solid phases
%--------------------------------------------------------------------------
figure(4)
clf
axis equal
vu_plot_n=vu_plot(:,mode)/norm(norm(vu_pf_Q(:,mode)))*a; % rescale
% vu_plot_n=sign*vu_plot(:,mode);
hold on
femplot(  x          ,   conns, 'Color', 'g');
femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
femplot(  x+vu_plot_n,   connss{2}, 'Color', 'r');
femplot(  x+vu_plot_n,   connss{3}, 'Color', 'b');
title("Mode shape of frequency "  + num2str(sfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------

elseif strcmp(design,'design4') || strcmp(design,'design8')

% plot the  displacement
%--------------------------------------------------------------------------
figure(4)
clf
sign=1; % control modes shape phase ( 0 or 180 degrees)
axis equal
hold on
vu_plot_n=sign*vu_plot(:,mode)/norm(norm(vu_pf_Q(:,mode)))*a; %rescale to see
% vu_plot_n=sign*vu_plot(:,mode);
femplot(  x,   conns, 'Color', 'g');
femplot(  x+vu_plot_n,   conns, 'Color', 'k');
title("Mode shape of frequency "  + num2str(sfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------

else

% plot the pressure solution and displacement
%--------------------------------------------------------------------------
figure(4)
clf
sign=1; % control modes shape phase ( 0 or 180 degrees)
scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),50,sign*vp(:,mode),'filled')
colorbar;
pmax = max(abs(vp_plot(:,mode)));
clim([-pmax pmax]);
colormap('jet');
axis equal
hold on
vu_plot_n=sign*vu_plot(:,mode)/norm(norm(vu_pf_Q(:,mode)))*a; %rescale to see
% vu_plot_n=sign*vu_plot(:,mode);
femplot(  x+vu_plot_n,   conns, 'Color', 'k');
title("Mode shape of frequency "  + num2str(sfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------

end

%% STEADY: PRESCRIBED DISPLACEMENT, PRESSURE AND ITS GRADIENTS

% reference position vector
%--------------------------------------------------------------------------
ones_p = ones(length(p) ,1);
xR=0*e1;
% ones_pp = ones(length(pp),1);
%--------------------------------------------------------------------------


% prescribed solid points are [corners_s]
%--------------------------------------------------------------------------
    % MODEL INPUT
    uM     = 0*(e1+e2);
    graduM = +0.0 *e1*e1 + 0.01   *e1*e2 +...
              0.01 *e2*e1 - 0.0   *e2*e2;
    % first-order comp. homogenization, prescribed disp. at the corners
    % note input index of x() is always "mesh" type.
    % delta xs_p
    dxs_p = x(p)-xR*ones_p;
    u_p    =  uM *ones_p + dot(graduM.', dxs_p);
    % transform to equivalent scalar (no need for sparsicity on a vector)
    su_p   = full(tm2sm_v(u_p, ee, 1));

%--------------------------------------------------------------------------




% stack both in a hybrid state vector containing disp. and pressure
wP = [su_p];

%--------------------------------------------------------------------------



% %% STATIONARY: SOLVE RIGHT/LEFT LINEAR SYSTEM OF EQUATIONS

%--------------------------------------------------------------------------
% right constrain modes
S_F_P = - inv(mu_F_F  ) * mu_F_P  ;
% right stationary hybrid state vector
su_f  =  S_F_P * wP;

% left constrain modes (Y appears always transposed along the expressions)
Y_F_P = - inv(mu_F_F.') * mu_P_F.';
% left stationary hybrid state vector
wl_f  =  Y_F_P * wP;
%--------------------------------------------------------------------------



% %% STATIONARY: ASSIGN DISP. SOLUTION TO TENSOR VECTOR
%--------------------------------------------------------------------------
% 'ured' is  reduced displacement thanks
% to the elimitation of the dependend points due to periodic bc
%
%  length(ured) = length([dofs_pc; dofs_un; dofs_in])
%  length(u   ) = length([dofs_c ; dofs_un; dofs_in; dofs_de])
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% free (reduced) displacements is the right stationary solution
    % declare
    % nodes are between index 1 to 2*length(f)
    % u_f = wr_f(1:2:2*length(f)-1)*e1+ wr_f(2:2:2*length(f))*e2;
     su_pf_Q = [su_p
                  su_f];

     vu_pf_Q = su_pf_Q(1:2:end)*e1+ su_pf_Q(2:2:end)*e2;

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% declare displacement tensor vectors (total # of dofs)
% assign to vu. Note vu is in [nodes_s] order
vu_pf_Q([p f])  = vu_pf_Q;
% free unused variables
clear vu_aux; clear ured;

% prepare for plot: length(vu_plot)>length(vu)
    % declare disp. vector for plot -  whole mesh
    vu_plot = zeros(1, ee,   n, 1);
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot(  nodes_s) = vu_pf_Q;
    % assign macro displacement to fluid nodes
%     vu_plot(nodes_f) = uM*ones(nf,1);
%--------------------------------------------------------------------------




disp('codeblock: PLOT')

% %% STATIONARY: PLOT PRESSURE/DISP SHAPE (RESCALED DISP.)


% plot the pressure solution and displacement
%--------------------------------------------------------------------------
if ns~=0 && nf~=0
    figure(5)
    clf
    % scatter(dot(  x,e1),dot(  x,e2),30,vp_plot,'filled')
    scatter(dot( x(nodes_f),e1),dot(x(nodes_f),e2),30,vp,'filled')
    colorbar;
    pmax = max(abs(vp_plot));
    caxis([-pmax pmax]);
    colormap('jet');
    axis equal
    hold on
    vu_plot_n= vu_plot/norm(norm(vu_pf_Q))*a; % rescale
    % vu_plot_n= vu_plot; % dont't rescale
    femplot(  x+vu_plot_n,   conns, 'Color', 'k');
    % title("Stationary problem (check if normalized)" );
    title("Stationary (check if normalized), x_R=(" ...
           +string(dot(xR,e1))+","+string(dot(xR,e2))+")" );
    hold off
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
if nf==0
    figure(5)
    clf
    axis equal
    hold on
    % vu_plot_n= vu_plot/norm(norm(vu)); % rescale
    vu_plot_n= vu_plot; % dont't rescale
    femplot(  x,   conns, 'Color', 'k');
    femplot(  x+vu_plot_n,   conns, 'Color', 'r');
    title("Stationary (check if normalized), x_R=(" ...
           +string(dot(xR,e1))+","+string(dot(xR,e2))+")" );
    hold off
end
%--------------------------------------------------------------------------




%% REDUCED COUPLED DYNAMIC MODEL
% Comments
%--------------------------------------------------------------------------
% 1. check about transposition .' or hermitian transposition '
% regular transposition is being used.
% 2. notation is, for instance, tlam_P_P means matrix with second order
% tensor components of size = (P,P).
% t - matrix with second order tensor components
% v - matrix with first order (vector) tensor components 
% s - matrix with scalar components
%--------------------------------------------------------------------------



% compute reduced matrices on P nodes with Q eigenmodes
%--------------------------------------------------------------------------
% quasistatic part
lam_lw  =  full( ...
                 lam_P_P   +  lam_P_F * S_F_P   +   Y_F_P.' * lam_F_P  ...
                + Y_F_P.'  *  lam_F_F * S_F_P ...
                                                                         );

mu_lw   =  full( ...
                 mu_P_P   +   mu_P_F * S_F_P   +   Y_F_P.' *  mu_F_P  ...
                + Y_F_P.' *   mu_F_F * S_F_P  ...
                                                                         );

lam_Q_P =   psi_F_Q.' * lam_F_P  +   psi_F_Q.' * lam_F_F * S_F_P  ;

% lam_P_Q =   lam_P_F   * phi_F_Q  +     Y_F_P.' * lam_F_F * phi_F_Q;

lam_P_Q =   lam_Q_P.';



% dynamics part
I_Q_Q    =    eye(n_modes);
LAM_Q_Q  =    psi_F_Q(:,1:n_modes)'*mu_F_F*phi_F_Q(:,1:n_modes);
O_Q_2    =    zeros(n_modes,2);
%--------------------------------------------------------------------------


% convert equivalent scalar matrix back to tensor matrix form
%--------------------------------------------------------------------------

% tlam_p_p   = sm2tm(lam_qs(1:dim*length(p),1:dim*length(p)),ee,2);
% tmu_p_p   = sm2tm(mu_qs(1:dim*length(p),1:dim*length(p)),ee,2);

tlam_lw   = sm2tm(lam_lw,ee,2);
tmu_lw   = sm2tm(mu_lw,ee,2);

% lam_P_Q hybrid matrix partitioning (for constitutive relations)
% vm_p_Q   = sm2tm(lam_P_Q(1:dim*length(p)    ,:),ee,1);
% vm_Q_p   = sm2tm(lam_Q_P(: , 1:dim*length(p)).',ee,1).';
vm_p_Q   = sm2tm(lam_P_Q,ee,1);
vm_Q_p   = sm2tm(lam_Q_P.',ee,1).';

%--------------------------------------------------------------------------


% compute homogenized macroscopic solid stress   
%--------------------------------------------------------------------------
%
% !! still need to left transpose hEs and hAs !!
% this can be done after
%
% note: h stands for 'homogenized' and s for 'stress'
hsA = 1/(a^2)*    l3transpose(dxs_p.' * tlam_lw  * ones_p  ,ee,3) ;
hsB = 1/(a^2)*     ltranspose(dxs_p.' * tlam_lw  * dxs_p)  ;
% hsC = 1/(a^2)*                dxs_p.' * vlam_p_pp * ones_pp ;
% hsD = 1/(a^2)*                dxs_p.' * vlam_p_pp * dxf_pp  ;
hsE = 1/(a^2)*    l3transpose(dxs_p.' * tmu_lw   * ones_p  ,ee,3) ;
hsF = 1/(a^2)*     ltranspose(dxs_p.' * tmu_lw   * dxs_p)  ;
% hsG = 1/(a^2)*                dxs_p.' * vmu_p_pp * ones_pp  ;
% hsH = 1/(a^2)*                dxs_p.' * vmu_p_pp * dxf_pp   ;
hsL = 1/(a*a)*               (vm_p_Q.'* dxs_p)              ;%corrected 20/04/2023-
% for the evolution equation
hsLs = 1/(a*a)*              (vm_Q_p  * dxs_p )             ;%corrected 20/04/2023-
hsl  =                       (vm_Q_p  * dxs_p )             ;

% convert to matrix
% tC_homog = tm2sm_v(hsF,ee,4);

% full contraction of each tensor for an estimation of order
abshsA = sqrt(  dddot(hsA,hsA) );
abshsB = sqrt( ddddot(hsB,hsB) );
% abshsC = sqrt(   ddot(hsC,hsC) );
% abshsD = sqrt(  dddot(hsD,hsD) );
abshsE = sqrt(  dddot(hsE,hsE) );
abshsF_BIOT = sqrt( ddddot(hsF,hsF) );
% abshsG_BIOT = sqrt(   ddot(hsG,hsG) );
% abshsH = sqrt(  dddot(hsH,hsH) );

abshsL=[];
for i=1:size(hsL,1)
abshsL =   [abshsL  sqrt( ddot(hsL(i),hsL(i)) )];
end
abshsL = abs(abshsL); 

abshsLs=[];
for i=1:size(hsLs,1)
abshsLs =   [abshsLs  sqrt( ddot(hsLs(i),hsLs(i)) )];
end
abshsLs = abs(abshsLs); 


% compute homogenized macroscopic solid momentum rate
%--------------------------------------------------------------------------
% still need to left transpose hEm and hAm
% note: h stands for 'homogenized' and m for 'momentum'
hmA = 1/(a^2)*                ones_p.' * tlam_lw  * ones_p  ;
hmB = 1/(a^2)*                ones_p.' * tlam_lw  * dxs_p   ;
% hmC = 1/(a^2)*                ones_p.' * vlam_p_pp * ones_pp ;
% hmD = 1/(a^2)*                ones_p.' * vlam_p_pp * dxf_pp  ;
hmE = 1/(a^2)*                ones_p.' * tmu_lw   * ones_p  ;
hmF = 1/(a^2)*                ones_p.' * tmu_lw   * dxs_p   ;
% hmG = 1/(a^2)*                ones_p.' * vmu_p_pp  * ones_pp ;
% hmH = 1/(a^2)*                ones_p.' * vmu_p_pp  * dxf_pp  ;
hmL = 1/(a*a)*               (ones_p.' * vm_p_Q).'           ;%corrected 20/04/2023-
% for the evolution equation
hmLs= 1/(a*a)*               (vm_Q_p  * ones_p)             ; %corrected 20/04/2023-
hml =                        (vm_Q_p  * ones_p)             ;


% full contraction of each tensor for an estimation of order
abshmA_BIOT = sqrt(    ddot(hmA,hmA) );
abshmB = sqrt(   dddot(hmB,hmB) );
% abshmC = sqrt(     dot(hmC,hmC) );
% abshmD = sqrt(    ddot(hmD,hmD) );
abshmE = sqrt(    ddot(hmE,hmE) );
abshmF = sqrt(   dddot(hmF,hmF) );
% abshmG = sqrt(     dot(hmG,hmG) );
% abshmH_BIOT = sqrt(    ddot(hmH,hmH) );

abshmL = [];
for i=1:size(hmL,1)
abshmL =  [abshmL  sqrt( dot(hmL(i),hmL(i)) )];
end

abshmLs = [];
for i=1:size(hmLs,1)
abshmLs =  [abshmLs  sqrt( dot(hmLs(i),hmLs(i)) )];
end


%% WEIGHTED DENSITIES

% %--------------------------------------------------------------------------
% % computed with analytical expression
% r3 = 5e-3;
% r2 = 7.5e-3;
% a  = 20e-3;
% V3 = pi*r3^2;
% V2 = pi*r2^2-V3;
% V1 = a^2-V2-V3;
% 
% rho_weighted_ava=(V1*matProp{1}.rho+V2*matProp{2}.rho+V3*matProp{3}.rho)/a^2;
% 
% 
% % computed with area of mesh
% rho_weighted_av=(V{1}*matProp{1}.rho+V{2}*matProp{2}.rho+V{3}*matProp{3}.rho)/a^2;
% %--------------------------------------------------------------------------




% %% PLOT COEFFICIENTS' MAGNITUDE VALUE
% clear etaticklabels
% 
% % build the x tick labels
% % mode labels
% etaticklabels = cell(1,n_modes);
% for i=1:n_modes
%     etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
%                          ' Hz)$ \eta_{' num2str(i) '} $'];
% end     
% 
% % make artificial x axis for plot
% xx=1:8+n_modes;
% 
% figure(6)
% clf
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% 
% % SOLID MOMENTUM
% %--------------------------------------------------------------------------
% subplot(2,1,1)
% Ym = [abshmA_BIOT abshmB abshmC abshmD abshmE abshmF ...
%       abshmG abshmH_BIOT abshmL];
% Xm = {'(BIOT) $ \ddot{\vec{u}}_M $', ...
%           ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
%           ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
%           ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
%           ' $ \vec{u}_{_{M}} $', ...
%           ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
%           ' $\mathrm{p}_{_{M}} $', ...
%           ' (BIOT) $\vec{\nabla} \mathrm{p}_{_{M}} $'};
% Xm = [Xm etaticklabels];
% plot(xx,log10(Ym), 'kx');
% hold on;
% plot(xx([1 8]),log10(Ym([1 8])), 'bx');
% plot(xx(9:8+n_modes),log10(Ym(9:8+n_modes)), 'gx');
% mincoeff = min(Ym(1),Ym(8)); % A and H are Biot coeffs
% plot(log10(mincoeff*ones(1,8)), "Color", "b");
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xm), 'XTickLabels',  Xm);
% title( '$\dot{\vec{\pi}}_{_{M}} $ Solid Momentum coefficients', 'interpreter', 'latex')
% ylabel('log10( absolute     value ) ')
% axi = gca;
% axi.FontSize = 16; 
% ylim([-10 4])
% %--------------------------------------------------------------------------
% clear Xm  mincoeff;
% 
% 
% % SOLID STRESS
% %--------------------------------------------------------------------------
% subplot(2,1,2)
% Ys = [abshsA abshsB abshsC abshsD abshsE abshsF_BIOT abshsG_BIOT ...
%       abshsH abshsL];
% Xs = {' $ \ddot{\vec{u}}_M $', ...
%       ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
%       ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $ \vec{u}_{_{M}} $', ...
%       ' (BIOT) $\vec{\nabla} \vec{u}_{_{M}} $', ...
%       ' (BIOT) $\mathrm{p}_{_{M}} $', ...
%       ' $\vec{\nabla} \mathrm{p}_{_{M}} $'};
% Xs = [Xs etaticklabels];
% plot(xx,log10(Ys), 'kx');
% hold on;
% mincoeff = min(Ys(6),Ys(7)); % F and G are Biot coeffs
% plot(xx([6 7]),log10(Ys([6 7])), 'bx');
% plot(xx(9:8+n_modes),log10(Ys(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), "Color", "b");
% title( '$\sigma_{_{M}} $ Solid Stress coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xs), 'XTickLabels',  Xs);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% ylim([-10 10])
% %--------------------------------------------------------------------------
% clear Xs  mincoeff;


% % FLUID VELOCITY
% %--------------------------------------------------------------------------
% subplot(2,2,3)
% Yv = [abshvA_BIOT abshvB abshvC abshvD abshvE abshvF abshvG ...
%       abshvH_BIOT abshvL];
% Xv = {'(BIOT) $ \ddot{\vec{u}}_M $', ...
%       ' $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
%       ' $ \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $ \vec{u}_{_{M}} $', ...
%       ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
%       ' $\mathrm{p}_{_{M}} $', ...
%       ' (BIOT) $\vec{\nabla} \mathrm{p}_{_{M}} $'};
% Xv = [Xv etaticklabels];
% plot(xx,log10(Yv), 'kx');
% hold on;
% mincoeff = min(Yv(1),Yv(8)); % A and H are Biot coeffs
% plot(xx([1 8]),log10(Yv([1 8])), 'bx');
% plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), "Color", "b");
% title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xv  mincoeff;
% 
% 
% 
% % FLUID VOLUMETRIC STRAIN
% %--------------------------------------------------------------------------
% subplot(2,2,4)
% Ye = [absheA absheB_BIOT absheC_BIOT absheD absheE absheF absheG ...
%       absheH absheL];
% Xe = {' $ \ddot{\vec{u}}_M $', ...
%       ' (BIOT) $\vec{\nabla} \ddot{\vec{u}}_{_{M}} $', ...
%       ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $\vec{\nabla} \ddot{\mathrm{p}}_{_{M}} $', ...
%       ' $ \vec{u}_{_{M}} $', ...
%       ' $\vec{\nabla} \vec{u}_{_{M}} $', ...
%       ' $\mathrm{p}_{_{M}} $', ...
%       ' $\vec{\nabla} \mathrm{p}_{_{M}} $'};
% Xe = [Xe etaticklabels];
% plot(xx,log10(Ye), 'kx');
% hold on;
% mincoeff = min(Ye(2),Ye(3)); % B and C are Biot coeffs
% plot(xx([2 3]),log10(Ye([2 3])), 'bx');
% plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), "Color", "b");
% % legend('coeff. abs', 'smallest biot');
% title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xe  mincoeff;

clear etaticklabels axi;

% %% DIMENSIONLESS DRIVING FIELDS (CONSTITUTIVE EQUATIONS)
% 
% % characteristic values
% %--------------------------------------------------------------------------
% % l0 = ay;
% % l0 = thi;
% 
% p0=0;
% 
% % max min
% f0 = 2000;
% % f0 = 100;
% 
% % homogenized quasistatic density
% hrhos = dot(e1,hmA,e1);
% 
% l0 = min( sqrt((matProp{1}.kappa + 4/3*matProp{1}.G)/hrhos)/f0 , ...
%           sqrt((                       matProp{1}.G)/hrhos)/f0 ) ;    %WAVELENGTH?
% u0 = 1e-6*l0;   % small strain
% %--------------------------------------------------------------------------
% 
% % macroscopic solid stress absolute value
% %--------------------------------------------------------------------------
% dimless_abshsA      = u0*f0^2     * abshsA      ;
% dimless_abshsB      = u0*f0^2/l0  * abshsB      ;
% dimless_abshsC      = p0*f0^2     * abshsC      ;
% dimless_abshsD      = p0*f0^2/l0  * abshsD      ;
% dimless_abshsE      = u0          * abshsE      ;
% dimless_abshsF_BIOT = u0/l0       * abshsF_BIOT ;
% dimless_abshsG_BIOT = p0          * abshsG_BIOT ;
% dimless_abshsH      = p0/l0       * abshsH      ;
% dimless_abshsL      = f0^2        * abshsL      ;
% %--------------------------------------------------------------------------
% 
% 
% % macroscopic solid momentum rate absolute value
% %--------------------------------------------------------------------------
% dimless_abshmA_BIOT = u0*f0^2     *abshmA_BIOT  ;
% dimless_abshmB      = u0*f0^2/l0  *abshmB       ;
% dimless_abshmC      = p0*f0^2     *abshmC       ;
% dimless_abshmD      = p0*f0^2/l0  *abshmD       ;
% dimless_abshmE      = u0          *abshmE       ;
% dimless_abshmF      = u0/l0       *abshmF       ;
% dimless_abshmG      = p0          *abshmG       ;
% dimless_abshmH_BIOT = p0/l0       *abshmH_BIOT  ;
% dimless_abshmL      = f0^2        *abshmL       ;
% %--------------------------------------------------------------------------
% 
% 
% % macroscopic fluid acceleration (v dot) absolute value
% %--------------------------------------------------------------------------
% dimless_abshvA_BIOT = u0*f0^2     *abshvA_BIOT  ;
% dimless_abshvB      = u0*f0^2/l0  *abshvB       ;
% dimless_abshvC      = p0*f0^2     *abshvC       ;
% dimless_abshvD      = p0*f0^2/l0  *abshvD       ;
% dimless_abshvE      = u0          *abshvE       ;
% dimless_abshvF      = u0/l0       *abshvF       ;
% dimless_abshvG      = p0          *abshvG       ;
% dimless_abshvH_BIOT = p0/l0       *abshvH_BIOT  ;
% dimless_abshvL      = f0^2        *abshvL       ;
% %--------------------------------------------------------------------------
% 
% % macroscopic fluid volumetric rate rate absolute value
% %--------------------------------------------------------------------------
% dimless_absheA      = u0*f0^2     *absheA       ;
% dimless_absheB_BIOT = u0*f0^2/l0  *absheB_BIOT  ;
% dimless_absheC_BIOT = p0*f0^2     *absheC_BIOT       ;
% dimless_absheD      = p0*f0^2/l0  *absheD       ;
% dimless_absheE      = u0          *absheE       ;
% dimless_absheF      = u0/l0       *absheF       ;
% dimless_absheG      = p0          *absheG  ;
% dimless_absheH      = p0/l0       *absheH       ;
% dimless_absheL      = f0^2        *absheL       ;
% %--------------------------------------------------------------------------
% 
% 
% 
% % PLOT DIMENSIONLESS COEFFICIENTS' MAGNITUDE VALUE
% 
% % mode labels
% etaticklabels = cell(1,n_modes);
% for i=1:n_modes
%     etaticklabels{i} =  ['(' num2str(round(freqs(i),0)) ...
%                          ' Hz)$ \eta_{' num2str(i) '}^* $'];
% end  
% 
% % make artificial x axis for plot
% xx=1:8+n_modes;
% 
% % without dimensions
% 
% figure(7)
% clf
% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% sgtitle([ 'DIMENSIONLESS :'...
%          '$ f_0 $ = ' num2str(f0) ...
%          '[Hz], $l_0$ =' num2str(l0*100,'%.1f') ...
%          '[cm], $ p_0 $ =' num2str(p0) '[Pa]'], 'Interpreter','latex');  
% % SOLID MOMENTUM
% %--------------------------------------------------------------------------
% subplot(2,1,1)
% Ym = [dimless_abshmA_BIOT dimless_abshmB dimless_abshmC ... 
%       dimless_abshmD dimless_abshmE dimless_abshmF ...
%       dimless_abshmG dimless_abshmH_BIOT dimless_abshmL];
% Xm = {'(BIOT) $ \ddot{\vec{u}}_M^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
%       ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $ \vec{u}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
%       ' $\mathrm{p}_{_{M}}^{*}  $', ...
%       ' (BIOT) $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% Xm = [Xm etaticklabels];
% plot(xx,log10(Ym), 'kx');
% hold on;
% mincoeff = min(Ym(1),Ym(8)); % A and H are Biot coeffs
% plot(xx([1 8]),log10(Ym([1 8])), 'bx');
% plot(xx(9:8+n_modes),log10(Ym(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), 'b');
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xm), 'XTickLabels',  Xm);
% title( '$\dot{\vec{\pi}}_{_{M}} $ Solid Momentum coefficients', 'interpreter', 'latex')
% ylabel('log10( absolute     value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xm Ym mincoeff;
% 
% 
% % SOLID STRESS
% %--------------------------------------------------------------------------
% subplot(2,1,2)
% Ys = [dimless_abshsA dimless_abshsB dimless_abshsC ... 
%       dimless_abshsD dimless_abshsE dimless_abshsF_BIOT ...
%       dimless_abshsG_BIOT dimless_abshsH dimless_abshsL];
% Xs = {' $ \ddot{\vec{u}}_M^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
%       ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
%       ' $ \vec{u}_{_{M}}^{*}  $', ...
%       ' (BIOT) $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
%       ' (BIOT) $\mathrm{p}_{_{M}}^{*}  $', ...
%       ' $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% Xs = [Xs etaticklabels];
% plot(xx,log10(Ys), 'kx');
% hold on;
% mincoeff = min(Ys(6),Ys(7)); % F and G are Biot coeffs
% plot(xx([6 7]),log10(Ys([6 7])), 'bx');
% plot(xx(9:8+n_modes),log10(Ys(9:8+n_modes)), 'gx');
% plot(log10(mincoeff*ones(1,8)), 'b');
% title( '$\sigma_{_{M}} $ Solid Stress coefficients', 'interpreter', 'latex')
% hold off;
% grid on;
% set(gca, 'XTick', 1:length( Xs), 'XTickLabels',  Xs);
% ylabel('log10( absolute value ) ')
% axi = gca;
% axi.FontSize = 16; 
% %--------------------------------------------------------------------------
% clear Xs Ys mincoeff;
% 
% 
% % % FLUID VELOCITY
% % %--------------------------------------------------------------------------
% % subplot(2,2,3)
% % Yv = [dimless_abshvA_BIOT dimless_abshvB dimless_abshvC ...
% %       dimless_abshvD dimless_abshvE dimless_abshvF dimless_abshvG ...
% %       dimless_abshvH_BIOT dimless_abshvL];
% % Xv = {'(BIOT) $ \ddot{\vec{u}}_M^{*}  $', ...
% %       ' $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
% %       ' $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
% %       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
% %       ' $ \vec{u}_{_{M}} ^{*} $', ...
% %       ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*} $', ...
% %       ' $\mathrm{p}_{_{M}}^{*} $', ...
% %       ' (BIOT) $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% % Xv = [Xv etaticklabels];
% % plot(xx,log10(Yv), 'kx');
% % hold on;
% % mincoeff = min(Yv(1),Yv(8)); % A and H are Biot coeffs
% % plot(xx([1 8]),log10(Yv([1 8])), 'bx');
% % plot(xx(9:8+n_modes),log10(Yv(9:8+n_modes)), 'gx');
% % plot(log10(mincoeff*ones(1,8)), 'b');
% % % legend('coeff. abs', 'smallest biot');
% % title( '$\dot{\vec{v}}_{_{M}} $ Fluid Velocity coefficients', 'interpreter', 'latex')
% % hold off;
% % grid on;
% % set(gca, 'XTick', 1:length( Xv), 'XTickLabels',  Xv);
% % ylabel('log10( absolute value ) ')
% % axi = gca;
% % axi.FontSize = 16; 
% % %--------------------------------------------------------------------------
% % clear Xv Yv mincoeff;
% % 
% % 
% % 
% % % FLUID VOLUMETRIC STRAIN
% % %--------------------------------------------------------------------------
% % subplot(2,2,4)
% % Ye = [dimless_absheA dimless_absheB_BIOT dimless_absheC_BIOT ...
% %       dimless_absheD dimless_absheE dimless_absheF ...
% %       dimless_absheG dimless_absheH dimless_absheL];
% % Xe = {' $ \ddot{\vec{u}}_M^{*}  $', ...
% %       ' (BIOT) $\vec{\nabla}^{*}  \ddot{\vec{u}}_{_{M}}^{*}  $', ...
% %       ' (BIOT) $ \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
% %       ' $\vec{\nabla}^{*}  \ddot{\mathrm{p}}_{_{M}}^{*}  $', ...
% %       ' $ \vec{u}_{_{M}}^{*}  $', ...
% %       ' $\vec{\nabla}^{*}  \vec{u}_{_{M}}^{*}  $', ...
% %       ' $\mathrm{p}_{_{M}}^{*}  $', ...
% %       ' $\vec{\nabla}^{*}  \mathrm{p}_{_{M}}^{*}  $'};
% % Xe = [Xe etaticklabels];
% % plot(xx,log10(Ye), 'kx');
% % hold on;
% % mincoeff = min(Ye(2),Ye(3)); % B and C are Biot coeffs
% % plot(xx([2 3]),log10(Ye([2 3])), 'bx');
% % plot(xx(9:8+n_modes),log10(Ye(9:8+n_modes)), 'gx');
% % plot(log10(mincoeff*ones(1,8)), 'b');
% % % legend('coeff. abs', 'smallest biot');
% % title( '$\ddot{e}_{_{M}} $ Fluid Volumetric Strain coefficients', 'interpreter', 'latex')
% % hold off;
% % grid on;
% % set(gca, 'XTick', 1:length( Xe), 'XTickLabels',  Xe);
% % ylabel('log10( absolute value ) ')
% % axi = gca;
% % axi.FontSize = 16; 
% % %--------------------------------------------------------------------------
% % clear Xe Ye mincoeff;
% 
% clear etaticklabels axi;


% PLOT SOLID STIFFNESS TENSOR

% C4  = tm2sm_v(hsF,ee,4);
% C4v = t2voigt(C4 ,ee,4)
% A1 = C4v(1,1)
% B1 = C4v(1,2)
% hE = (A1-B1)*(A1+2*B1)/(A1+B1)
% hnu= B1/(A1+B1)
% PRINT HOMOGENIZED MATERIAL PARAMETERS

% % ISOTROPIC artificially created (-> get the first component)
% 
% % mode resonance
% moderes=2;
% 
% % local resonance material parameters in constitutive relations
% vas = (a^2)* hmL(moderes)
% tbs = (a^2)* hsL(moderes)
% % scs = (a^2)* heL(moderes)
% % vds = (a^2)* hvL(moderes)
% 
% % local resonance material parameters in evolution equation
% valphas = (a^2)* hmLs(moderes)
% tbetas  = (a^2)* hsLs(moderes)
% % sgammas = (a^2)* heLs(moderes)
% % vdeltas = (a^2)* hvLs(moderes)

%% DISPERSION BLOCH ANALYSIS (COMSOL)

% figure(21)
scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),20,'k*');
hold on
grid on
% ylim([0 2000])
% colormap("cool");
% colorbar('Ticks',[0.01,0.99],...
%          'TickLabels',{'v', 'u'})
% caxis([0 1])

% % plot
% %--------------------------------------------------------------------------
% figure(24)
% hold on
% grid on
% for i=1:number_of_waves
% scatter(range_of_dimless_k,imag(frequencies(i,:)),'r.')
% end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% hold off
% %--------------------------------------------------------------------------

%% DISPERSION CURVE REGULAR SOLID TERMS (QUASISTATIC)
disp('DISPERSION CALCULATION')


isorthogonal = [];

% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    % Density matrix (k)
        % first line of M
        Muu=tm2sm_v(hmA,ee,2);
 
       
    % Elasticity matrix (k)
      
        Kuu=tm2sm_v(dot(k,hsF,k),ee,2);
                                                           
    
    % density matrix 
    M=Muu;
    % elasticity matrix 
    K=Kuu;
    
    % compute eigevalues
    
    % lambda(:,i)=sort(eig(K,M));    



      % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------


    % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(eigenvectors(:,1)'*eigenvectors(:,2))  ];
    %--------------------------------------------------------------------------
    




    i=i+1;
    end
%--------------------------------------------------------------------------


% % plot
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% % figure(23)
% hold on
% % grid on
% for i=1:number_of_waves
% % scatter(k_plot,real(frequencies(i,:)),'g.') %warning('for gammaXM')
% scatter(range_of_dimless_k,real(frequencies(i,:)),'g.')
% % scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
% end
% % title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% % ylim([0 2000])
% % hold off
% %--------------------------------------------------------------------------


% plot eigenvector orthogonal?
%--------------------------------------------------------------------------
figure
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
end
colorbar
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


% plot log10 imag part of lambda
%--------------------------------------------------------------------------
figure
hold on
grid on
log10_im_lambda = log10(abs(imag(lambda)));
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'b.')
scatter(range_of_dimless_k,log10_im_lambda(i,:), '.')
end
%--------------------------------------------------------------------------

%% DISPERSION CURVE REGULAR SOLID TERMS (QUASISTATIC) + m-MODES   


isorthogonal = [];


% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % % wavenumber
    % if(dimless_k<0)
    % k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    % elseif(dimless_k<1)
    % k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    % else
    % k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    % end
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end
    
    % Density matrix (k)
    % first line of M
    Muu=tm2sm_v(hmA,ee,2);
    Mue=tm2sm_v(hmL.',ee,1);
    
    
    % third line of M
    Meu=tm2sm_v(hml.',ee,1).' ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix (k)
    % first line of K
    Kuu=tm2sm_v(dot(k,hsF,k),ee,2);
    Kue= O_Q_2.' ;
    
    % third line of K
    Keu= O_Q_2 ;
    Kee= LAM_Q_Q ;                                                          

% density matrix 
M=[Muu Mue;
   Meu Mee];
% elasticity matrix 
K=[Kuu Kue;
   Keu Kee];
    
    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------
     
    % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(eigenvectors(:,1)'*eigenvectors(:,2) + ...
                         eigenvectors(:,1)'*eigenvectors(:,3) + ...
                         eigenvectors(:,1)'*eigenvectors(:,4) + ...
                         eigenvectors(:,1)'*eigenvectors(:,5) + ...
                         eigenvectors(:,2)'*eigenvectors(:,3) + ...
                         eigenvectors(:,2)'*eigenvectors(:,4) + ...
                         eigenvectors(:,2)'*eigenvectors(:,5) + ...
                         eigenvectors(:,3)'*eigenvectors(:,4) + ...
                         eigenvectors(:,3)'*eigenvectors(:,5) + ...
                         eigenvectors(:,4)'*eigenvectors(:,5)) ];
    %--------------------------------------------------------------------------

    % % are the eigenvectors orthonormal?
    % %--------------------------------------------------------------------------
    % isorthogonal = [isorthogonal ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,2))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,3))'*real(eigenvectors(:,4))) + ...
    %                  abs(real(eigenvectors(:,3))'*real(eigenvectors(:,5))) + ...
    %                  abs(real(eigenvectors(:,4))'*real(eigenvectors(:,5))) ];
    % %--------------------------------------------------------------------------

    
    i=i+1;
    end
%--------------------------------------------------------------------------


% plot one color
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
% hold on
% grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'b.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'g.')
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


% % plot eigenvector orthogonal?
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% % figure(23)
% hold on
% grid on
% for i=1:number_of_waves
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
% end
% colorbar
% % title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% % ylim([0 2000])
% % hold off
% %--------------------------------------------------------------------------


% % plot log10 imag part of lambda
% %--------------------------------------------------------------------------
% % figure
% hold on
% grid on
% log10_im_lambda = log10(abs(imag(lambda)));
% for i=1:number_of_waves
% % scatter(k_plot,real(frequencies(i,:)),'b.')
% scatter(range_of_dimless_k,log10_im_lambda(i,:), '.')
% end
% ylim([0 2000])
% %--------------------------------------------------------------------------


%% DISPERSION CURVE REGULAR SOLID TERMS (QUASISTATIC) +  s-MODES 
% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    % Density matrix 
    % first line of M
    Muu=tm2sm_v(hmA,ee,2);
    Mue=      +1i*tm2sm_v(dot(k,hsL.'),ee,1);
 
    % third line of M
    Meu=      +1i*tm2sm_v(dot( - hsl.',k),ee,1).'  ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(dot(k,hsF,k),ee,2);
    Kue= O_Q_2.' ;
    
    % third line of K
    Keu= O_Q_2 ;
    Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=[Muu Mue;
       Meu Mee];
    % elasticity matrix 
    K=[Kuu Kue;
       Keu Kee];
    
    lambda(:,i)=sort(eig(K,M));       
    
    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
% hold on
% grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'k.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'k.')
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------




%%  DISPERSION CURVE REGULAR SOLID TERMS (QUASISTATIC) + m-MODES  + s-MODES 
% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    % Density matrix 
    % first line of M
    Muu=tm2sm_v(hmA,ee,2);
    Mue=tm2sm_v(hmL.',ee,1)           +1i*tm2sm_v(dot(k,hsL.'),ee,1);
 
    % third line of M
    Meu=tm2sm_v(+ hml.',ee,1).'      +1i*tm2sm_v(dot( - hsl.',k),ee,1).'  ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(dot(k,hsF,k),ee,2);
    Kue= O_Q_2.' ;
    
    % third line of K
    Keu= O_Q_2 ;
    Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=[Muu Mue;
       Meu Mee];
    % elasticity matrix 
    K=[Kuu Kue;
       Keu Kee];
    
    lambda(:,i)=sort(eig(K,M));       
    
    i=i+1;
    end
%--------------------------------------------------------------------------


% plot
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
% hold on
% grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'k.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'k.')
end
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


%% DISPERSION CURVE ALL QUASISTATIC TERMS

isorthogonal = [];

% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    % Density matrix 
    % first line of M
    Muu=tm2sm_v(hmA+dot(k,hsB,k),ee,2)+1i*tm2sm_v(dot(-hmB,k)+dot(k,hsA),ee,2);
    % Mue=tm2sm_v(hmL.',ee,1)           +1i*tm2sm_v(dot(k,hsL.'),ee,1);
 
    % third line of M
    % Meu=tm2sm_v(+ hml.',ee,1).'      +1i*tm2sm_v(dot( - hsl.',k),ee,1).'  ;
    % Mee= I_Q_Q                                                            ;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(hmE+dot(k,hsF,k),ee,2)+1i*tm2sm_v(dot(-hmF,k)+dot(k,hsE),ee,2);
    % Kue= O_Q_2.' ;
    
    % third line of K
    Keu= O_Q_2 ;
    Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=Muu;
    % elasticity matrix 
    K=Kuu;
    
       % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------


    % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(eigenvectors(:,1)'*eigenvectors(:,2))  ];
    %--------------------------------------------------------------------------
    % % are the eigenvectors orthonormal?
    % %--------------------------------------------------------------------------
    % isorthogonal = [isorthogonal ...
    %                  abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2)))  ];
    % %--------------------------------------------------------------------------      
    
    i=i+1;
    end
%-----------------------------------------------------------------------


% % plot
% %--------------------------------------------------------------------------
% figure
% frequencies = sqrt(lambda)/2/pi;
% c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% % figure(23)
% hold on
% % grid on
% for i=1:number_of_waves
% % scatter(k_plot,real(frequencies(i,:)),'b.')
% scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
% end
% % title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% ylim([0 2000])
% hold off
% %--------------------------------------------------------------------------

% plot eigenvector orthogonal?
%--------------------------------------------------------------------------
figure
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
end
colorbar
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


% plot log10 imag part of lambda
%--------------------------------------------------------------------------
figure
hold on
grid on
log10_im_lambda = log10(abs(imag(lambda)));
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'b.')
scatter(range_of_dimless_k,log10_im_lambda(i,:), '.')
end
%--------------------------------------------------------------------------

%% DISPERSION CURVE ALL TERMS


% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------

isorthogonal = [];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.005;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    % Density matrix 
    % first line of M
    Muu=tm2sm_v(hmA+dot(k,hsB,k),ee,2)+1i*tm2sm_v(dot(-hmB,k)+dot(k,hsA),ee,2);
    Mue=tm2sm_v(hmL.',ee,1)           +1i*tm2sm_v(dot(k,hsL.'),ee,1);
 
    % third line of M
    Meu=tm2sm_v(+ hml.',ee,1).'      +1i*tm2sm_v(dot( - hsl.',k),ee,1).'  ;
    Mee= I_Q_Q                                                            ;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(hmE+dot(k,hsF,k),ee,2)+1i*tm2sm_v(dot(-hmF,k)+dot(k,hsE),ee,2);
    Kue= O_Q_2.' ;
    
    % third line of K
    Keu= O_Q_2 ;
    Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=[Muu Mue;
       Meu Mee];
    % elasticity matrix 
    K=[Kuu Kue;
       Keu Kee];
    
    % lambda(:,i)=sort(eig(K,M));       
    
    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------

    
    % % are the eigenvectors orthonormal?
    % %--------------------------------------------------------------------------
    % isorthogonal = [isorthogonal ...
    %                  abs(eigenvectors(:,1)'*eigenvectors(:,2) + ...
    %                      eigenvectors(:,1)'*eigenvectors(:,3) + ...
    %                      eigenvectors(:,1)'*eigenvectors(:,4) + ...
    %                      eigenvectors(:,1)'*eigenvectors(:,5) + ...
    %                      eigenvectors(:,2)'*eigenvectors(:,3) + ...
    %                      eigenvectors(:,2)'*eigenvectors(:,4) + ...
    %                      eigenvectors(:,2)'*eigenvectors(:,5) + ...
    %                      eigenvectors(:,3)'*eigenvectors(:,4) + ...
    %                      eigenvectors(:,3)'*eigenvectors(:,5) + ...
    %                      eigenvectors(:,4)'*eigenvectors(:,5)) ];
    % %--------------------------------------------------------------------------
     
     % are the eigenvectors orthonormal?
    %--------------------------------------------------------------------------
    isorthogonal = [isorthogonal ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,2))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,1))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,3))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,2))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,3))'*real(eigenvectors(:,4))) + ...
                     abs(real(eigenvectors(:,3))'*real(eigenvectors(:,5))) + ...
                     abs(real(eigenvectors(:,4))'*real(eigenvectors(:,5))) ];
    %--------------------------------------------------------------------------



    i=i+1;
    end
%-----------------------------------------------------------------------

% 
% plot one color only
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'r.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'g.')
end
colorbar
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------


% % plot eigenvector orthogonal?
% %--------------------------------------------------------------------------
% frequencies = sqrt(lambda)/2/pi;
% c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% % figure(23)
% hold on
% grid on
% for i=1:number_of_waves
% scatter(range_of_dimless_k,real(frequencies(i,:)),6,log10(isorthogonal),'filled' );
% end
% colorbar
% % title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
% % ylim([0 2000])
% % hold off
% %--------------------------------------------------------------------------

 %% DISPERSION CURVE SPACE HARMONIC HOMOGENIZATION

% initialize vector with phases due to Bloch boundary condition
e_p = zeros(length(p),1);

% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------

isorthogonal = [];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------




% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    

    % assign phases for given k
    %--------------------------------------------------------------------------
    e_p = exp(1i*dot(dxs_p,k)); % 
    %--------------------------------------------------------------------------
    
    % homogenized coefficients for given k
    %--------------------------------------------------------------------------
    C_pi = 1/a^2*   e_p'* tlam_lw*e_p;
    D_pi = 1/a^2*   e_p'* tmu_lw *e_p;
    b    = 1/a^2*  (e_p'* vm_p_Q).';
    %--------------------------------------------------------------------------

    % mass matrix 
    Muu=tm2sm_v(C_pi,ee,2);  Mue= tm2sm_v(b',ee,1);
    Meu=(a^2)*Mue'       ;  Mee= I_Q_Q;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(D_pi,ee,2);  Kue= O_Q_2.' ;
    Keu= O_Q_2 ;             Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=[Muu Mue;
       Meu Mee];
    % elasticity matrix 
    K=[Kuu Kue;
       Keu Kee];
    
    % lambda(:,i)=sort(eig(K,M));       
    
    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------

 


    i=i+1;
    end
%-----------------------------------------------------------------------

% 
% plot one color only
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'r.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
scatter(range_of_dimless_k,imag(frequencies(i,:)),'b.')
end
colorbar
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------

 %% DISPERSION CURVE BLOCH MODE SYNTHESIS


T = [ phi_F_Q S_F_P; zeros(lenght(p),lenght(f)) ; eye(lenght(p))  ];

% initialize vector with phases due to Bloch boundary condition
e_p = zeros(length(p),1);

% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------

isorthogonal = [];

% number of depend. variables describing the continuum
% (u1,u2,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 2+n_modes; % propagating 
%--------------------------------------------------------------------------

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.01;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = 0:kstep:1;
%--------------------------------------------------------------------------

% wave number normalized by pi/ax
%--------------------------------------------------------------------------
Gamma_index=   round(length(range_of_dimless_k)/3)+1;
X_index    = 2*round(length(range_of_dimless_k)/3)+1;

% dimless wavenumber normalized by pi/ax
k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
               range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
               range_of_dimless_k(X_index+1:end    )    *ay/ax ];
%--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------




% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    

    % assign phases for given k
    %--------------------------------------------------------------------------
    e_p = exp(1i*dot(dxs_p,k)); % 
    %--------------------------------------------------------------------------
    
    % homogenized coefficients for given k
    %--------------------------------------------------------------------------
    C_pi = 1/a^2*   e_p'* tlam_lw*e_p;
    D_pi = 1/a^2*   e_p'* tmu_lw *e_p;
    b    = 1/a^2*  (e_p'* vm_p_Q).';
    %--------------------------------------------------------------------------

    % mass matrix 
    Muu=tm2sm_v(C_pi,ee,2);  Mue= tm2sm_v(b',ee,1);
    Meu=(a^2)*Mue'       ;  Mee= I_Q_Q;

    % Elasticity matrix
    % first line of K
    Kuu=tm2sm_v(D_pi,ee,2);  Kue= O_Q_2.' ;
    Keu= O_Q_2 ;             Kee= LAM_Q_Q ;                                                                  

    
    % density matrix 
    M=[Muu Mue;
       Meu Mee];
    % elasticity matrix 
    K=[Kuu Kue;
       Keu Kee];
    
    % lambda(:,i)=sort(eig(K,M));       
    
    % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % sort d and V with respect to eigenvalues
    %--------------------------------------------------------------------------
    [omega2,Ind] = sort(diag(eigenvalues));
    eigenvectors    = eigenvectors(:,Ind);
    %--------------------------------------------------------------------------

    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = omega2;
    %--------------------------------------------------------------------------

 


    i=i+1;
    end
%-----------------------------------------------------------------------

% 
% plot one color only
%--------------------------------------------------------------------------
frequencies = sqrt(lambda)/2/pi;
c_x = 2*pi*(frequencies(:,Gamma_index+1)-frequencies(:,Gamma_index)) / (kstep*pi/ax); % speed of sound gamma-X
% figure(23)
hold on
grid on
for i=1:number_of_waves
% scatter(k_plot,real(frequencies(i,:)),'r.')
scatter(range_of_dimless_k,real(frequencies(i,:)),'r.')
scatter(range_of_dimless_k,imag(frequencies(i,:)),'b.')
end
colorbar
% title("slope \Gamma-X is c = " + num2str(c_x)+"[m/s]")
ylim([0 2000])
% hold off
%--------------------------------------------------------------------------