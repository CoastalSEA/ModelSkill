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
    usertools = {'Add Properties','Delete Properties','Tabulate Properties',...
                 'Plot Properties','Plot Hypsometry','Plot X-Y property sets'};
   [selection,ok] = listdlg('Liststring',usertools,...
                                 'PromptString','Select tool:',...  
                                 'ListSize',[180,100],...
                                 'SelectionMode','single');
    if ok==0, return; end  %user cancelled selection
    src.Text = usertools{selection};
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    
    
    switch src.Text
        case 'Add Properties'
            GD_ImportData.gridMenuOptions(mobj,src,gridclasses);
%             addFormProperties(obj,muicat); 
        case 'Delete Properties'
            GD_ImportData.gridMenuOptions(mobj,src,gridclasses);
%             delFormProperties(obj,muicat);      
        case 'Tabulate Properties'
            getGrossPropsTable(mobj,true);
        case 'Plot Properties'
            getPropsPlot(mobj);
        case 'Plot Hypsometry'
            getHypsPlot(mobj);
        case 'Plot X-Y property sets'
            getXYPlot(mobj)
    end
end
%%
function tabledata = getGrossPropsTable(mobj,istab)
    %tabulate the gross properties in a figure
%     rN = 15; %No. of rows in the gross properties struct
    muicat = mobj.Cases;
    promptxt = 'Select a Model:';
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
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    hold on
    select = 1; count = 1;
    while select>0
        [tabledata,varnum,ok] = selectPlotData(mobj);
        if ok~=0
            data = tabledata.DataTable{:,varnum};
            x = tabledata.RowNames;            
            legtext = inputdlg('Legend text','Plot properties',1,{num2str(count)});
            plot(ax,x,data,'LineWidth',1,'DisplayName',legtext{1});
            if count==1
                xlabel('Time step (years)')
                ylabel(tabledata.VariableLabels{varnum});
            end
            
            answer = questdlg('Add another dataset?','Plot properties','Yes');
            if strcmp(answer,'No')
                select = 0;
                legend
            else
                count = count+1;
            end
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
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    hold on
    select = 1; count = 1;
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to plot:','Select timestep:'};
    while select>0
        [cobj,~,irow] = selectCaseDatasetRow(mobj.Cases,[],...
                                                 gridclasses,promptxt,1);
        if ~isempty(cobj)
            wl = cobj.Data.WaterLevels.DataTable(irow,:);
            hyps = cobj.Data.Hypsometry.DataTable(irow,:);  
            zcentre = cobj.Data.Hypsometry.Dimensions.Z;
            zcentre(end,1) = zcentre(end,1)+0.1; %crude offset to highest point so that HW is visible
            zsurf = hyps.SurfaceArea;
            zvol = hyps.Volume;
            casedesc = cobj.Data.Form.Description;
            labs = sprintf('%s Area at %s',casedesc,hyps.Properties.RowNames{1});
            plot(zsurf,zcentre,'--','DisplayName',labs);
            labv = sprintf('%s Volume at %s',casedesc,hyps.Properties.RowNames{1});
            plot(zvol,zcentre,'-.','DisplayName',labv);

            if count==1                
                %add water levels at mouth
                lightgrey = [0.7,0.7,0.7];
                h1 = plot(xlim, wl.zhw(1)*[1 1],':','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                h1 = plot(xlim, wl.zmt(1)*[1 1],'--','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
                h1 = plot(xlim, wl.zlw(1)*[1 1],':','Color',lightgrey);
                h1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
            end
        else 
            select = 0;
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
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    hold on
    select = 1; count = 1;
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
            end
            
            if count==1
                xlabel(tabledataX.VariableLabels{varnumX});
                ylabel(tabledataY.VariableLabels{varnumY});
            end
            %
            if strcmp(xname,yname)
                metatxt = xname;
            else
                metatxt = springtf('%s-%s',xname,yname);
            end
            
            plot(ax,dataX,dataY,'LineWidth',1,'DisplayName',metatxt);
            xpts = tabledataX.RowNames;    
            if ~isempty(xpts)
                txtlabel = ['\leftarrow ',char(string(xpts(1)))];
                text(dataX(1),dataY(1),txtlabel,'FontSize',8);
                txtlabel = ['\leftarrow ',char(string(xpts(end)))];
                text(dataX(end),dataY(end),txtlabel,'FontSize',8);
            end
            answer = questdlg('Add another dataset?','Plot properties','Yes');
            if strcmp(answer,'No')
                select = 0;
                legend
            else
                count = count+1;
            end
        else
            select = 0;
        end
    end
    hold off    
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

     
                 


