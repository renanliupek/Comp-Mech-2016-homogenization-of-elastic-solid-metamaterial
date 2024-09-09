function mesh=read_mphtxt(s)

%%%% FIRST VERSION, UNUSED %%%%%% 


% tensor basis defintion
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(b{1}, b{2});
e1 = ee(1);
e2 = ee(2);



%% cell x
%--------------------------------------------------------------------------
% define cell
mesh.x{1,1}='coordinate';
mesh.x{1,2}='tensor vec';

% import coordinates   
        % Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 6);
        
        % Specify range and delimiter
        opts.DataLines = [22, 10505];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "Created", "Var3", "Var4", "Var5", "Var6"];
        opts.SelectedVariableNames = ["VarName1", "Created"];
        opts.VariableTypes = ["double", "double", "string", "string", "string", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");
        
        % Import the data
        x_design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);
    

        % assign cell content
        mesh.x{2,1} = table2array(x_design);
        mesh.x{2,2} = mesh.x{2,1}(:,1)*e1+mesh.x{2,1}(:,2)*e2;
        
        % Clear temporary variables
        clear opts
%--------------------------------------------------------------------------

%% cell edge
%--------------------------------------------------------------------------
mesh.edge{1,1} = 'conn';
mesh.edge{1,2} = 'tag';

% import conn
        % Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 6);
        
        % Specify range and delimiter
        opts.DataLines = [10561, 11338];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "Created", "Var3", "Var4", "Var5", "Var6"];
        opts.SelectedVariableNames = ["VarName1", "Created"];
        opts.VariableTypes = ["double", "double", "string", "string", "string", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");
        
         % Import the data
        edge_design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

        % Convert to output type
        mesh.edge{2,1} = table2array(edge_design);

        % Clear temporary variables
        clear opts


% import tag
        % Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 6);
        
        % Specify range and delimiter
        opts.DataLines = [11342, 12119];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6"];
        opts.SelectedVariableNames = "VarName1";
        opts.VariableTypes = ["double", "string", "string", "string", "string", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");
        
        % Import the data
        etag_design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

        % Convert to output type
        mesh.edge{2,2} = table2array(etag_design);
        
        % Clear temporary variables
        clear opts

%--------------------------------------------------------------------------


%% cell elem
%--------------------------------------------------------------------------
mesh.elem{1,1} = 'conn';
mesh.elem{1,2} = 'tag';

% import conn
        % Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 6);
        
        % Specify range and delimiter
        opts.DataLines = [12129, 22409];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "Created", "by", "COMSOL", "Var5", "Var6"];
        opts.SelectedVariableNames = ["VarName1", "Created", "by", "COMSOL"];
        opts.VariableTypes = ["double", "double", "double", "double", "string", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, ["Var5", "Var6"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var5", "Var6"], "EmptyFieldRule", "auto");
        opts = setvaropts(opts, ["by", "COMSOL"], "TrimNonNumeric", true);
        opts = setvaropts(opts, ["by", "COMSOL"], "ThousandsSeparator", ",");
        
        % Import the data
        conn_design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

        % Convert to output type
        mesh.elem{2,1} = table2array(conn_design);
        
        % Clear temporary variables
        clear opts


% import tag
        % Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 6);
        
        % Specify range and delimiter
        opts.DataLines = [22413, Inf];
        opts.Delimiter = " ";
        
        % Specify column names and types
        opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6"];
        opts.SelectedVariableNames = "VarName1";
        opts.VariableTypes = ["double", "string", "string", "string", "string", "string"];
        
        % Specify file level properties
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        
        % Specify variable properties
        opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");
        
        % Import the data
        tag_design = readtable(['C:\GitLab\Tools\coho\meshes\' s '.txt'], opts);

        % Convert to output type
        mesh.elem{2,2} = table2array(tag_design);
        
        % Clear temporary variables
        clear opts

%--------------------------------------------------------------------------






end