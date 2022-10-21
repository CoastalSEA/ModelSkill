function getInletTools(mobj)
%
%-------function help------------------------------------------------------
% NAME
%   getInletTools.m
% PURPOSE
%   tools to extract and plot morphological properties of inlets and
%   channels
% USAGE
%   getInletTools(obj,classrec,muicat)
% INPUTS
%   mobj - ModelUI instance
% OUTPUT
%   Add Properties: add form properties to a gridded data set.
%   Delete Properties: delete ALL property tables associated with a  selected gridded data set.
%   Tabulate Properties: tabulate the gross properties for all timesteps in Selected Case.
%   Plot Properties: plot selected gross properties for selected Cases as a function of time.
%   Plot Hypsometry: plot surface area and volume hypsometry for one or more selected Cases and time intervals.
%   Plot X-Y property sets: plot selected gross properties as X and Y to show co-variant change 
% NOTES
%   uses tools in GDinterface 
% SEE ALSO
%   getUserTools. Called in ModelSkill
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2022
%--------------------------------------------------------------------------
%
    ok = 1;
    while ok>0
        usertools = {'Add Properties','Delete Properties',...
                     'Edit Inlet Definition','Plot Case Properties',...
                     'Tabulate Set Properties, f(t)','Plot Set Properties, f(t)',...
                     'Plot Hypsometries','Plot X-Y property sets, f(t)',...
                     'Plot rates of change, f(t)',...
                     'Plot Thalwegs','User Function'};
       [selection,ok] = listdlg('Liststring',usertools,...
                                     'PromptString','Select tool (Cancel to Quit):',...  
                                     'ListSize',[180,160],...
                                     'SelectionMode','single');
        if ok==0, continue; end  %user cancelled selection
        src.Text = usertools{selection};
        gridclasses = {'GD_ImportData'}; %add other classes if needed

        switch src.Text
            case 'Add Properties'
                FGDinterface.formMenuOptions(mobj,src,gridclasses);
            case 'Delete Properties'
                FGDinterface.formMenuOptions(mobj,src,gridclasses);
            case 'Plot Case Properties'
                getCasePropsFigure(mobj);
                ok = 0; %quit usertools selection to allow access to figure UI
            case 'Edit Inlet Definition'
                FGDinterface.formMenuOptions(mobj,src,{'GD_ImportData'});
            case 'Tabulate Set Properties, f(t)'
                getGrossPropsTable(mobj,true);
            case 'Plot Set Properties, f(t)'
                getPropsPlot(mobj);
            case 'Plot Hypsometries'
                getHypsPlot(mobj);
            case 'Plot X-Y property sets, f(t)'
                getXYPlot(mobj);
            case 'Plot rates of change, f(t)'
                getRatesPlot(mobj);
            case 'Plot Thalwegs'
                getThalwegs(mobj);
            case 'User Function'
                ms_userfunction(mobj);
        end
    end
end
%%
function getCasePropsFigure(mobj)
    %create figure for cf_section plot output as used in ChannelForm model
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to plot:','Select timestep:'};
    [cobj,~,irec] = selectCaseDatasetRow(mobj.Cases,[],...
        gridclasses,promptxt,1);
    if isempty(cobj) || isempty(irec), return; end
    %generate table and plot for display on Properties tab
    if ~isfield(cobj.Data,'GrossProps') || isempty(cobj.Data.GrossProps)
        getdialog('No Form Properties available for selected grid')
        return;
    end
    T = getDSTable(cobj.Data.GrossProps,irec,[]);
    
    hf = figure('Name','Gross properties', ...
        'Units','normalized', ...
        'Tag','PlotFig');
    hf.MenuBar = 'none';   %remove menu bar
    hf.ToolBar = 'figure'; %allow access to data tips and save tools

    %generate table of gross properties
    uitable('Parent',hf,'Data',T.DataTable{:,:},...
        'ColumnName',T.VariableNames,...
        'ColumnWidth',{80},...
        'RowName',T.DataTable.Properties.RowNames,...
        'Units','Normalized','Position',[0.05,0.84,0.9,0.14],...
        'Tag','grossprops');
    %user popup to select a type of plot
    popup = findobj(hf,'Style','popup','Tag','PlotList');
    if isempty(popup)
        plotlist = {'Hypsommetry','Cross-sections','Thalweg + Plan width',...
            'Along-channel properties','Hydraulic depth',...
            'Prism-Area ratio','Elevation-Area histogram',...
            'a/h and Vs/Vc','Hydraulics','Transgression'};
        uicontrol('Parent',hf,'Style','text',...
            'Units','Normalized','Position', [0.05 0.79 0.1 0.04],...
            'String', 'Select plot:');
        popup = uicontrol('Parent',hf,'Style','popup',...
            'String',plotlist,'Tag','PlotList',...
            'Units','Normalized','Position', [0.05 0.74 0.9 0.05],...
            'Callback', @(src,evdat)gd_property_plots(cobj,irec,src)); %#ok<NASGU>

        %Create push button to copy data to clipboard
        uicontrol('Parent',hf,'Style','pushbutton',...
            'String','>Table','UserData',T,'Tag','CopyButton',...
            'TooltipString','Copy grossproperties table to clipboard',...
            'Units','normalized','Position',[0.75 0.793 0.10 0.044],...
            'Callback',@copydata2clip);

        %create push button to create tab plot as a stand alone figure
        uicontrol('Parent',hf,'Style','pushbutton',...
            'String','>Figure','Tag','FigButton',...
            'TooltipString','Create plot as stand alone figure',...
            'Units','normalized','Position',[0.86 0.793 0.10 0.044],...
            'Callback',@(src,evdat)gd_property_plots(cobj,irec,src));

    else
        %update obj in Callbacks and table in UserData
        popup.Callback = @(src,evdat)gd_property_plots(cobj,irec,src);
        hb = findobj(hf,'Tag','CopyButton');
        hb.UserData = T;
        hf = findobj(hf,'Tag','FigButton');
        hf.Callback = @(src,evdat)gd_property_plots(cobj,irec,src);
        gd_property_plots(cobj,irec,popup); %set plot to current popup selection
    end
end
%%
function tabledata = getGrossPropsTable(mobj,istab)
    %tabulate the gross properties in a figure
    muicat = mobj.Cases;
    promptxt = 'Select a Model (Cancel to Quit):';
    [caserec,ok] = selectCase(muicat,promptxt,'single',0);     
    if ok<1, tabledata = []; return; end
    
    %retrieve additional data from data set
    cobj = getCase(mobj.Cases,caserec);
    if isfield(cobj.Data,'GrossProps')
        tabledata = cobj.Data.GrossProps;
    else
        warndlg('No GrossProps table for selected case');
        tabledata = []; 
        return;
    end
    
    figtitle = 'Gross properties';
    casedesc = mobj.Cases.Catalogue.CaseDescription(caserec);
    casedesc = sprintf('%s using %s',figtitle,casedesc);    
    if istab
        tablefigure(figtitle,casedesc,tabledata);
    end
end    
%%
function getPropsPlot(mobj)
    %plot gross properties from groups of results
    hf = figure('Name','Gross properties', ...
                'Units','normalized', ...
                'Tag','PlotFig');
    ax = axes(hf);
    lines = {'-','--',':','-.'};
    
    hold on
    select = 1; count = 0;
    while select>0
        [tabledata,varnum,ok] = selectPlotData(mobj);
        if ok~=0
            data = tabledata.DataTable{:,varnum};
            x = tabledata.RowNames;            
            legtext = inputdlg('Legend text','Plot properties',1,{num2str(count)});
            plot(ax,x,data,'LineStyle',lines{rem(count,4)+1},...
                                'LineWidth',1,'DisplayName',legtext{1});
            if count==0
                xlabel('Time step (years)')
                ylabel(tabledata.VariableLabels{varnum});
            end
            
            [select,count] = addAnotherDataset(count);
        else
            select = 0;
        end
    end
    hold off    
end
%%
function getHypsPlot(mobj)
    %plot the grid hypsometry
    hf = figure('Name','Gross properties', ...
                'Units','normalized', ...
                'Tag','PlotFig');
    ax = axes(hf);
    hold on
    select = 1; count = 1;
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to plot:','Select timestep:'};
    while select>0
        [cobj,~,irow] = selectCaseDatasetRow(mobj.Cases,[],...
                                                 gridclasses,promptxt,1);
        if isempty(cobj) || isempty(irow)
            select = 0;
        else
            wl = cobj.Data.WaterLevels.DataTable(irow,:);
            hyps = cobj.Data.Hypsometry.DataTable(irow,:);  
            zcentre = cobj.Data.Hypsometry.Dimensions.Z;
            zcentre(end,1) = zcentre(end,1)+0.1; %crude offset to highest point so that HW is visible
            zsurf = hyps.SurfaceArea;
            zvol = hyps.Volume;
            casedesc = cobj.Data.Grid.Description;
            labs = sprintf('%s Area at %s',casedesc,hyps.Properties.RowNames{1});
            plot(ax,zsurf,zcentre,'--','DisplayName',labs);
            labv = sprintf('%s Volume at %s',casedesc,hyps.Properties.RowNames{1});
            plot(ax,zvol,zcentre,'-.','DisplayName',labv);

            if count==1                
                %add water levels at mouth
                lightgrey = [0.7,0.7,0.7];
                h1 = plot(ax,xlim, wl.zhw(1)*[1 1],':','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                h1 = plot(ax,xlim, wl.zmt(1)*[1 1],'--','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                h1 = plot(ax,xlim, wl.zlw(1)*[1 1],':','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            end  
        end
    end
    hold off

    xlabel('Volume (m^3) and Area (m^2)'); 
    ylabel('Elevation (mAD)');
    legend('Location','SouthEast');
    title('System hypsometry','FontWeight','normal','FontSize',10);
end
%%
function getXYPlot(mobj)
    %plot gross properties from groups of results
    hf = figure('Name','Gross properties', ...
                'Units','normalized', ...
                'Tag','PlotFig');
    ax = axes(hf);
    lines = {'-','--',':','-.'};
    
    hold on
    select = 1; count = 0;
    while select>0
        [tabledataX,varnumX,okX] = selectPlotData(mobj);
        [tabledataY,varnumY,okY] = selectPlotData(mobj);
        if okX~=0 && okY~=0
            dataX = tabledataX.DataTable{:,varnumX};                    
            xname = tabledataX.Description;
            dataY = tabledataY.DataTable{:,varnumY};        
            yname = tabledataY.Description;
            if height(tabledataX)~=height(tabledataY)
                warndlg('Selected data are not the same length')
                select = 0;
                continue;
            end
            
            if count==0
                xlabel(tabledataX.VariableLabels{varnumX});
                ylabel(tabledataY.VariableLabels{varnumY});
            end
            %
            if strcmp(xname,yname)
                metatxt = xname;
            else
                metatxt = sprintf('%s-%s',xname,yname);
            end
            
            plot(ax,dataX,dataY,'LineStyle',lines{rem(count,4)+1},...
                       'Marker','.','LineWidth',1,'DisplayName',metatxt);
            xpts = tabledataX.RowNames;    
            if ~isempty(xpts)
                txtlabel = ['\leftarrow ',char(string(xpts(1)))];
                text(dataX(1),dataY(1),txtlabel,'FontSize',8);
                txtlabel = ['\leftarrow ',char(string(xpts(end)))];
                text(dataX(end),dataY(end),txtlabel,'FontSize',8);
            end
            
            [select,count] = addAnotherDataset(count);
        else
            select = 0;
        end
    end
    hold off    
end
%%
function getRatesPlot(mobj)
    %plot the rates of change of a selected variable
    hf = figure('Name','Gross properties', ...
                'Units','normalized', ...
                'Tag','PlotFig');
    ax = axes(hf);
    lines = {'-','--',':','-.'};
    
    hold on
    select = 1; count = 0;
    while select>0
        [tabledata,varnum,ok] = selectPlotData(mobj);
        if ok~=0
            vardata = tabledata.DataTable{:,varnum};  
            t = tabledata.RowNames; 
            varate = diff(vardata)./years(diff(t));
            varate = [0;varate];          
            legtext = inputdlg('Legend text','Plot properties',1,{num2str(count)});
            plot(ax,t,varate,'LineStyle',lines{rem(count,4)+1},...
                                'LineWidth',1,'DisplayName',legtext{1});
            if count==0
                xlabel('Time step (years)')
                varname = tabledata.VariableDescriptions{varnum};
                ylabel(sprintf('Rate of change of %s',varname));
            end
            
            [select,count] = addAnotherDataset(count);
        else 
            select = 0;
        end
    end
    hold off
end
%%
function getThalwegs(mobj)
    %plot the thalwegs (deepest point in channel) between user defined
    %start and end points
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt1 = {'Select Case to plot (Cancel to quit):','Select timestep:'};
    [obj,~,irec] = selectCaseDatasetRow(mobj.Cases,[],...
                                                 gridclasses,promptxt1,1);
    if isempty(obj) || isempty(irec), return; end
    desc = sprintf('%s at %s',obj.Data.Grid.Description,char(obj.Data.Grid.RowNames(irec)));
    grid = getGrid(obj,irec);   %grid for selected year
    [X,Y] = meshgrid(grid.x,grid.y);
    N = numel(X);
    xy = [reshape(X,[N,1]),reshape(Y,[N,1])];
    Z = grid.z';
    
    %get maximum water level to define 
    promptxt2 = {'Maximum accessible water level?','Depth exponent'};
    defaults = {num2str(max(Z,[],'all')),'2'};
    answer = inputdlg(promptxt2,'Water level',1,defaults);
    if isempty(answer), return; end %user cancelled
    maxwl = str2double(answer{1});
    dexp = str2double(answer{2});

    %accessible map (water) and use -Z as the cost map
    water = true(size(Z));
    water(isnan(Z) | Z>maxwl) = false;
 
    %cline = gd_selectpoints(grid,true);     
    promptxt3 = {'Select start of path','Select end of path'};
    gridmasked = grid;        gridmasked.z(~water') = NaN;
    points = gd_selectpoints(gridmasked,2,promptxt3,true);
    if any(isnan([points(:).x])), return; end
    
    %index of nearest grid point to selected start end end points    
    start = dsearchn(xy,[points(1).x,points(1).y]); 
    goal = dsearchn(xy,[points(2).x,points(2).y]);
    
    %find the shortest path taking account of the cost (depths)
    %Z(Z>maxwl) = 0;
    costfnc = @(z) -(min(z,[],'all')-z).^dexp; %weighted inverse depth to favour staying deep
    thalweg = a_star(water, costfnc(Z), start, goal);
    [idy,idx] = ind2sub(size(Z),thalweg);
    
    %plot base map of initial grid selection and defined mask
    hf = figure('Name','Thalwegs','Units','normalized','Tag','PlotFig');                            
    ax = gd_plotgrid(hf,grid);
    colormap(ax,'gray');
    lines = {'-','--',':','-.'};
    
    hs = findobj(ax.Children,'Type','surface');
    hs.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
    hp = plot(ax,[points(:).x],[points(:).y],'ok','MarkerSize',8,'MarkerFaceColor','w','Tag','mypoints');
    hp.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    plot(ax,grid.x(idx),grid.y(idy),'r','LineStyle',lines{1},'LineWidth',1,...
                                                      'DisplayName',desc);
    hold off
    title('Thalwegs between defined start and end points')
    legend
    
    select = 1; count = 1;
    while select>0
        [obj,~,irec] = selectCaseDatasetRow(mobj.Cases,[],...
                                                 gridclasses,promptxt1,1);                                            
        if isempty(obj) || isempty(irec)
            select = 0;   %user cancelled
        else
            desc = sprintf('%s at %s',obj.Data.Grid.Description,char(obj.Data.Grid.RowNames(irec))); 
            grid = getGrid(obj,irec);   %grid for selected year
            Z = grid.z';
            %update water mask to reflect new grid
            water = true(size(Z));
            water(isnan(Z) | Z>maxwl) = false;            
            %find the shortest path taking account of the cost (depths)
            %Z(Z>maxwl) = 0;
            thalweg = a_star(water, costfnc(Z), start, goal);
            if isempty(thalweg)
                hw = warndlg(sprintf('No path found for %s',desc));
                uiwait(hw)
                continue;
            end
            [idy,idx] = ind2sub(size(Z),thalweg);
            hold on
            plot(ax,grid.x(idx),grid.y(idy),'r','LineStyle',...
                 lines{rem(count,4)+1},'LineWidth',1,'DisplayName',desc);
            hold off 
            count = count+1;
        end                                       
    end
end
%%
function [tabledata,varnum,ok] = selectPlotData(mobj)
    %select data set and variable to be used
    tabledata = getGrossPropsTable(mobj,false); 
    if isempty(tabledata)
        varnum = []; ok = 0;
        return; 
    end
    PromptText = 'Select Variable'; 
    ListSize = [300,250];           
    [varnum,ok] = listdlg('Name','Variable', ...
            'ListSize',ListSize,...
            'PromptString',PromptText, ...
            'SelectionMode','single', ...
            'ListString',tabledata.VariableNames);  
end
%%
function [select,count] = addAnotherDataset(count)
    answer = questdlg('Add another dataset?','Plot properties','Yes','No','Yes');
    if strcmp(answer,'No')
        select = 0;
        legend
    else
        select = 1;
        count = count+1;
    end
end

     
                 


