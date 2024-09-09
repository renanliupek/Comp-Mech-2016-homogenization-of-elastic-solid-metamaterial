function mesh = read_mphtxt_54(s)

%--------------------------------------------------------------------------
% Renan Liupekevicius TU/e 30-05-2022
% READ_MPHTXT_2  import .txt file mesh to my matlab mesh
%   mesh = READ_MPHTXT_2(s)
%   
%   input  s    - string, file name
%   output mesh - mesh struct 
% 
% Assumptions
%  - volume element tag=1 for solid phase (to convert tag=3 to tag=1);
%  - Use of basis 'e1' and 'e2';
%  - Use COMSOL 5.6;
%  - 2D mesh;
%--------------------------------------------------------------------------


%% IMPORT NUMERIC PART

    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 9);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = " ";
    
    % Specify column names and types
    opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";

    % Specify variable properties
    opts = setvaropts(opts, ["VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"], "ThousandsSeparator", ",");

    
    % Import the data
    design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

    % Convert to output type
    num = table2array(design);
    
    % Clear temporary variables
    clear opts, clear design;

%% IMPORT STRING PART
% the string version is used to find the text delimiter within the mesh
% file separating coordinates from connectivity matrix etc.

   % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 1);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = "";
    
    % Specify column names and types
    opts.VariableNames = "CreatedByCOMSOLMultiphysics";
    opts.VariableTypes = "string";
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    
    % Specify variable properties
    opts = setvaropts(opts, "CreatedByCOMSOLMultiphysics", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "CreatedByCOMSOLMultiphysics", "EmptyFieldRule", "auto");
    
    % Import the data
    text = readmatrix(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

    % Clear temporary variables
    clear opts



%% SELECT COLUMN DELIMITERS

% the following imports the 2D quadrilateral mesh from COMSOL 5.4

% coordinate
    % find
    [iia,~] = find(text=='# Mesh vertex coordinates');
    if isempty(iia)
        [iia,~] = find(text=='# Mesh point coordinates');
    end

    % start index
    iia = iia + 1; % 1 lines below
    
    % end index
    iia(2)=iia + num(iia-4,1) - 1; % 4 lines above

% edge connectivity
    %find
    [iib,~] = find(text=='3 # number of vertices per element');
    
    % start index
    iib = iib + 3; % 3 lines below

    % end index
    iib(2) = iib + num(iib-2,1) - 1; % 2 lines above

% edge tag
    % start index
    iic    = iib(2) + 4; % 4 lines below

    % end index

    iic(2) = iic + num(iic-2,1) - 1; % 2 lines above

% elem connectivity
    % find
    [iid,~] = find(text=='9 # number of vertices per element');
    
    % start index
    iid = iid + 3; % 3 lines below

    % end index
    iid(2) = iid + num(iid-2) -1; % 2 lines above

% elem tag
    % start index
    iie = iid(2) + 4; % 4 lines below

    % end index
    iie(2) = iie + num(iie-2,1) - 1; % 2 lines above


%% DEFINE DATA STRUCTURE

% tensor basis defintion
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(b{1}, b{2});
e1 = ee(1);
e2 = ee(2);
clear b;

% x
mesh.x{1,1}='coordinate';
mesh.x{2,1}=num(iia(1):iia(2), [1 2]);
mesh.x{1,2}='tensor vec';
mesh.x{2,2}= mesh.x{2,1}(:,1)*e1 + mesh.x{2,1}(:,2)*e2;

% edge
mesh.edge{1,1} = 'conn';
mesh.edge{2,1} =  num(iib(1):iib(2), [1 2 3]);
mesh.edge{2,1} =  mesh.edge{2,1} + 1;
mesh.edge{1,2} = 'tag';
mesh.edge{2,2} =  num(iic(1):iic(2), [ 1 ]  );
mesh.edge{2,2} =  mesh.edge{2,2} + 1;

% elem
mesh.elem{1,1} = 'conn';
mesh.elem{2,1} =  num(iid(1):iid(2), [1 3 4 2 6 9 8 5 7] );
mesh.elem{2,1} =  mesh.elem{2,1} + 1;
mesh.elem{1,2} = 'tag';
mesh.elem{2,2} =  num(iie(1):iie(2), [ 1 ] );


%% DELETE CENTER NODES
% 
%     
%     % extract node number of a center element
%     icenter = mesh.elem{2,1}(:,9);
% 
%     % delete from tensor vector and matrix
%     mesh.x{2,1}(icenter,:) = [];
% 
    % delete from conn matrix
    mesh.elem{2,1}(:,9)    = [];
% 
% 
% % x
% mesh.x{1,2}='tensor vec';
% mesh.x{2,2}= mesh.x{2,1}(:,1)*e1 + mesh.x{2,1}(:,2)*e2;


end