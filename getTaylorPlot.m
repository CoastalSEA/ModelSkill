function obj = getTaylorPlot(mobj)
%
%-------function help------------------------------------------------------
% NAME
%   getTaylorPlot.m
% PURPOSE
%   call the Taylor plot function with either gridded or
%   timeseries data
% USAGE
%   obj = getTaylorPlot(mobj)
% INPUTS
%   mobj - handle to ModelUI class
% OUTPUT
%   obj - dummy value for use in DataManip
% NOTES
%   Taylor, K, 2001, Summarizing multiple aspects of model performance 
%   in a single diagram, JGR-Atmospheres, V106, D7.   
% SEE ALSO
%   ModelSkill.m, taylor_plot.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%
    obj = []; %dummy value for use in DataManip    
    gridrecs = find(strcmp(mobj.Cases.Catalogue.CaseClass,'MS_GridData'));
    tsrecs = find(strcmp(mobj.Cases.Catalogue.CaseClass,'muiUserData'));
    if ~isempty(gridrecs) && ~isempty(tsrecs)
        answer = questdlg('Select data type','Taylor','Grid','Timeseries','Grid');       
    elseif  ~isempty(gridrecs)
        answer = 'Grid';
    elseif ~isempty(tsrecs)
        answer = 'Timeseries';  
    else
        warndlg('No grid or timeseries data available')
        return
    end

    switch answer
        case 'Grid'
            get_Grid_TaylorPlot(mobj);
        case 'Timeseries'
            get_TS_TaylorPlot(mobj);
    end
end
%%
function get_Grid_TaylorPlot(mobj)
    %select timeseries data and call the taylor_plot function   
    %prompt user to select reference grid
    promptxt = 'Select a Reference grid:';
    [grid0,metatxt(1),ok] = grid_selection(mobj,promptxt);
    if ok<1, return; end
    
    %extract Skill parameters from MS_RunParameters class
    [skill,ok] = getSkillParameters(mobj,grid0.x,grid0.y);
    if ok<1, return; end

    %iteratively select test data set
    rLim = 2;
    answer = 'New';    
    while ~strcmp(answer,'Exit')
        promptxt = 'Select a Test grid';
        [grid,metatxt(2),ok] = grid_selection(mobj,promptxt);        
        if ok<1, return; end

        taylor_plot(grid0.z,grid.z,metatxt,answer,rLim,skill)

        answer = questdlg('Add or Delete data set?','Taylor plot',...
            'Add','Delete','Exit','Add');
        if strcmp(answer,'Add')
            addoption = questdlg('Use the same Reference Grid?','Taylor plot',...
                                 'Yes','No','Yes');
            if strcmp(addoption,'No')
                promptxt = 'Select a Reference grid:';
                [grid0,metatxt(1),ok] = grid_selection(mobj,promptxt);
                if ok<1, return; end
                [skill,ok] = getSkillParameters(mobj,grid0.x,grid0.y);
                if ok<1, return; end
            end
        end
    end
end
%%
function [grid,casedesc,ok] = grid_selection(mobj,promptxt)
    %prompt user to select grid
    grid = []; casedesc = {''};
    
    [caserec,ok] = selectRecord(mobj.Cases,'PromptText',promptxt,...
                         'CaseClass',{'MS_GridData'},'ListSize',[200,200]);                                                                                   
    if ok<1, return; end

    [grid,timetxt,caserec,~] = MS_GridData.getCaseGridData(mobj.Cases,caserec);
    %retrieve description from Results Case
    casedesc = mobj.Cases.Catalogue.CaseDescription(caserec);
    casedesc = {sprintf('%s at %s',casedesc,timetxt)};
end
%%
function [skill,ok] = getSkillParameters(mobj,x,y)
    %extract Skill parameters from MS_RunParameters class
    ok = 1; skill = [];    
    if isfield(mobj.Inputs,'MS_RunParams')
        runparams = mobj.Inputs.MS_RunParams;
        skill.Ro = runparams.maxcorr;
        skill.n  = runparams.skillexponent;
        skill.Inc = true;                    %flag to include skill score
        skill.W = runparams.skillwindow;
        subdomain = runparams.skillsubdomain;
        skill.SD = getSubDomain(x,y,subdomain);
        skill.iter = runparams.skilliteration;
    else
        warndlg('Run parameters for skill score have not been defined');
        ok = 0;
        return
    end
end
%%
function sd = getSubDomain(x,y,subdomain)
    %find the subdomain in integer grid indices defined by x,y range
    %subdomain defined as [x0,xN,y0,yN];
    
    if isempty(subdomain) || length(subdomain)~=4
        subdomain = [min(x),max(x),min(y),max(y)];
    end
    ix0 = find(x<=subdomain(1),1,'last');
    ixN = find(x>=subdomain(2),1,'first');
    iy0 = find(y<=subdomain(3),1,'last');
    iyN = find(y>=subdomain(4),1,'first');
    sd.x = [ix0,ix0,ixN,ixN];
    sd.y = [iyN,iy0,iy0,iyN];
end
%%
function get_TS_TaylorPlot(mobj)
    %select timeseries data adn call the taylor_plot function
    %prompt user to select reference timeseries
    promptxt = 'Select a Timeseries to use as Reference case?';
    [z0,metatxt(1),ok] = ts_selection(mobj,promptxt);
    if ok<1, return; end
    
    %extract Skill parameters from MS_RunParameters class
    [skill,ok] = getSkillParameters(mobj,z0,[]);
    if ok<1, return; end

    %iteratively select test data set
    rLim = 2;
    answer = 'New';    
    while ~strcmp(answer,'Exit')
        promptxt = 'Select a Timeseries to use as Test case?';
        [zi,metatxt(2),ok] = ts_selection(mobj,promptxt);        
        if ok<1, return; end

        taylor_plot(z0,zi,metatxt,answer,rLim,skill)

        answer = questdlg('Add or Delete data set?','Taylor plot',...
            'Add','Delete','Exit','Add');
        if strcmp(answer,'Add')
            addoption = questdlg('Use the same Reference Timeseries?','Taylor plot',...
                                 'Yes','No','Yes');
            if strcmp(addoption,'No')
                promptxt = 'Select a Timeseries to use as Reference case?';
                [z0,metatxt(1),ok] = ts_selection(mobj,promptxt);
                if ok<1, return; end
            end
        end        
    end  
end
%%
function [ts,metatxt,ok] = ts_selection(mobj,promptxt)
    %prompt user to select timeseries
    ts = []; metatxt = {''};
    [caserec,ok] = selectRecord(mobj.Cases,'PromptText',promptxt,...
                         'CaseClass',{'muiUserData'},'ListSize',[200,200]);                                                                                   
    if ok<1, return; end
    dst = getDataset(mobj.Cases,caserec,1);
    dst_names = dst.VariableNames;
    idx = 1;
    if length(dst_names)>1
        [idx,ok] = listdlg('Name','Variables', ...
                    'PromptString','Select a variable:', ...
                    'ListSize',[200,200],...
                    'SelectionMode','single', ...
                    'ListString',dst_names);
        if ok<1, return; end
    end
    metatxt = mobj.Cases.Catalogue.CaseDescription(caserec);
    ts = dst2tsc(dst,[],idx);
end

            
            
