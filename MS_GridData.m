classdef MS_GridData
%
%-------class help---------------------------------------------------------
% NAME
%   MS_GridData.m
% PURPOSE
%   Class to import grid data, adding the results to dstable
%   and a record in a dscatlogue (as a property of muiCatalogue)
% USAGE
%   obj = MS_GridData()
% SEE ALSO
%   uses dstable and dscatalogue
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%    
    properties  
        Data
    end
    
    properties (Hidden, SetAccess=private)
        CaseIndex       %case index assigned when class instance is loaded     
    end
    
    methods 
        function obj = MS_GridData()               
            %class constructor                          
        end
    end
%%    
    methods (Static)
        function loadData(muicat,~)
            %read and load a data set from a file
            % mobj - handle to modelui instance 
            % classname - name of class being loaded
            obj = MS_GridData; %initialise instance of class
            
            [fname,path,nfiles] = getfiles('MultiSelect','on',...
                'FileType','*.txt;*.grd;*.csv','PromptText','Select file(s):');
            if ~iscell(fname)
                fname = {fname};   %single select returns char
            end
            
            if nfiles>1
                timesteps = {num2str(0:1:nfiles-1)};
            else
                timesteps = {'0'};
            end
            timesteps = inputdlg('Define timesteps:','Load grid',1,timesteps);
            if isempty(timesteps), return; end
            timesteps = years(str2double(split(timesteps)));
            
            dsp = setDSproperties(obj);  %initialise dsproperties for data
            
            isflip = false;
            %check whether cartesian grid needs to be flipped 
            answer = questdlg('Reverse origin','Data import','No','Yes','No');
            if strcmp(answer,'Yes'), isflip = true; end

            %now load each file selected
            nhead = 1; %number of header lines in file   
            for i=1:nfiles
               [data,~,filename{i,1}] = readInputData(obj,[path,fname{i}],nhead);   %#ok<AGROW>
               if isempty(data), continue; end 
               griddata(i,:,:) = formatGridData(obj,data,isflip); %#ok<AGROW>
            end

            %load the results into a dstable  
            dst = dstable(griddata,'RowNames',timesteps,'DSproperties',dsp); 
            dst.Dimensions.X = unique(data(:,1));
            dst.Dimensions.Y = unique(data(:,2));

            %assign metadata about data
            dst.Source = filename;
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setDataSetRecord(obj,muicat,dst,'data');           
        end 
%%
        function [grid,timetxt,caserec,irec] = getCaseGridData(muicat,caserec)
            %get the selected grid
            if isempty(caserec)
                promptxt = 'Select a Model:';
                [caserec,~] = selectCase(muicat,promptxt,'single',0);     
                if isempty(caserec), grid = []; return; end
            end
            dst = getDataset(muicat,caserec,1);
            if height(dst)>1
                %propmpt user to select timestep
                list = dst.DataTable.Properties.RowNames;
                irec = listdlg('PromptString','Select timestep:',...
                               'Name','Tab plot','SelectionMode','single',...
                               'ListSize',[200,200],'ListString',list);
            else
                irec = 1;
            end

            grid.x = dst.Dimensions.X;
            grid.y = dst.Dimensions.Y;
            grid.z = squeeze(dst.Z(irec,:,:));
            timetxt = dst.DataTable.Properties.RowNames{irec};
        end
%%
        function reGridData(mobj,~,~)
            %regrid an imported data set to match another grid or to user
            %specified dimensions. Because the grid changes size need to
            %apply to all time steps and save as a new record
            muicat = mobj.Cases;
            %prompt user to select an existing form model to modify
            promptxt = 'Select a Model Grid to Regrid?';
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                        'CaseClass',{'MS_GridData'},'ListSize',[200,200]);
            if ok<1, return; end
            cobj = getCase(muicat,caserec);
            dst = cobj.Data.Dataset;
            x = dst.Dimensions.X;
            y = dst.Dimensions.Y;
            casedesc = muicat.Catalogue.CaseDescription(caserec);
            
            %prompt user to define source of X-Y axes
            answer = questdlg('Use Run Parameter Grid dimensions, or an existing Cartesian Grid?',...
                              'Regrid','New','Existing','New');                          
            if strcmp(answer,'New')
                %use Run Parameter definitions to 
                if isfield(mobj.Inputs,'MS_RunParams')
                    runparams = mobj.Inputs.MS_RunParams;                    
                    [xaxis,yaxis] = getplotaxes(cobj,runparams,x,y);
                    [X,Y] = meshgrid(xaxis,yaxis);  
                    xi = X(1,:)';
                    yi = Y(:,1); 
                    metatxt = sprintf('Regridded from %s using run parameter definition',casedesc);
                else                    
                    warndlg('Set Grid Dimensions in order to create a New grid');
                    return
                end
            else
                %select an existing Cartesian grid to use
                promptxt = 'Select a Model Grid to define grid X and Y ranges?';
                [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                         'CaseClass',{'MS_GridData'},'ListSize',[200,200]);                                                                                   
                if ok<1, return; end
                refdst = getDataset(muicat,caserec,1);
                griddesc = muicat.Catalogue.CaseDescription(caserec);
                metatxt = sprintf('Regridded from %s using %s',casedesc,griddesc);
                refgrid.x = refdst.Dimensions.X;
                refgrid.y = refdst.Dimensions.Y;
                [X,Y] = meshgrid(refgrid.x,refgrid.y);
                xi = X(1,:)';
                yi = Y(:,1);
            end   
            %
            Nx = length(xi); Ny = length(yi);
            for i=1:height(dst)
                z =squeeze(dst.DataTable.Z(i,:,:));
                Z = griddata(x,y,z',X,Y); %#ok<GRIDD>                
                impdata(i,:,:) = reshape(Z',1,Nx,Ny);                 %#ok<AGROW>
            end
            newdst = copy(dst);
            newdst.DataTable.Z = impdata;
            newdst.Dimensions.X = xi;
            newdst.Dimensions.Y = yi;

            %assign metadata about data
            newdst.Source = {metatxt};
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setDataSetRecord(cobj,muicat,newdst,'data');     
            mobj.DrawMap;
        end    
%%
        function subGridData(mobj,~,~)
            %interactively define a subgrid and save grid as a new case
            muicat = mobj.Cases;
            %prompt user to select an existing form model to modify
            promptxt = 'Select a Model Grid to Regrid?';
            [caserec,ok] = selectRecord(muicat,'PromptText',promptxt,...
                        'CaseClass',{'MS_GridData'},'ListSize',[200,200]);
            if ok<1, return; end
            obj = getCase(muicat,caserec);
            dst = obj.Data.Dataset;
            x = dst.Dimensions.X;
            y = dst.Dimensions.Y;
            z = squeeze(dst.Z(1,:,:)); %first row/time
            casedesc = dst.Description;
            

            %create subgrid selection plot
            [subdomain,defaults] = getSubDomain(obj,x,y,z);
            %extract selected grid
            [xo,yo,~,ixo,iyo] = getsubgrid(x,y,z',subdomain);
            Nx = length(xo); Ny = length(yo);
            for i=1:height(dst)
                zi =squeeze(dst.DataTable.Z(i,:,:))';
                zo = zi(min(iyo):max(iyo),min(ixo):max(ixo));
                impdata(i,:,:) = reshape(zo',1,Nx,Ny); %#ok<AGROW>
            end

            %extract defined subgrid for each row in table
            newdst = copy(dst);
            newdst.DataTable.Z = impdata;
            newdst.Dimensions.X = xo;
            newdst.Dimensions.Y = yo;

            %assign metadata about data
            newdst.Source = sprintf('Subgrid of %s using subdomain: %s',...
                casedesc,defaults);
            %setDataRecord classobj, muiCatalogue obj, dataset, classtype
            setDataSetRecord(obj,muicat,newdst,'data');
            mobj.DrawMap;
        end
    end   
%%
    methods
        function addData(obj,classrec,~,muicat) 
            %add additional data to an existing user dataset
            [fname,path,~] = getfiles('MultiSelect','off',...
                'FileType','*.txt;*.csv','PromptText','Select file:');
            nhead = 1;
            [data,~,~] = readInputData(obj,[path,fname],nhead);
            if isempty(data), return; end
            
            isflip = false;
            %check whether cartesian grid needs to be flipped 
            answer = questdlg('Reverse origin','Data import','No','Yes','No');
            if strcmp(answer,'Yes'), isflip = true; end
            
            dst = obj.Data.Dataset;      %selected dstable
            existimes = dst.RowNames;
            ok = 1;
            while ok>0
                timestep = inputdlg('Define timestep:','Load grid',1,{'X'});
                if isempty(timestep), return; end
                timestep = years(str2double(timestep));
                if ~ismember(existimes,timestep), ok=0; end %ensure value entered is unique
            end
            
            %reformat data and add as a new row to existing table
            griddata = formatGridData(obj,data,isflip);            
            nrec = height(dst);
            dst.DataTable{nrec+1,1} = griddata;
            dst.RowNames(nrec+1) = timestep;

            %assign metadata about data
            dst.Source{nrec+1} = fname;
            
            obj.Data.Dataset = dst;  
            updateCase(muicat,obj,classrec);
        end        
%%
        function deleteData(obj,classrec,catrec,muicat)
            %delete variable or rows from a dataset
            dst = obj.Data.Dataset;      %selected dstable
            
            delist = dst.VariableNames;
            %select variable to use
            promptxt = {sprintf('Select Variable')}; 
            att2use = 1;
            if length(delist)>1
                [att2use,ok] = listdlg('PromptString',promptxt,...
                                 'Name',title,'SelectionMode','single',...
                                 'ListSize',[250,100],'ListString',delist);
                if ok<1, return; end  
            end
            promptxt = sprintf('Delete variable: %s?',delist{att2use});
            selopt = questdlg(promptxt,'Delete variable',...
                                      'Yes','No','No');
            if strcmp(selopt,'Yes')
                dst.(delist{att2use}) = [];  %delete selected variable

                obj.Data.Dataset = dst;            
                updateCase(muicat,obj,classrec);
                getdialog(sprintf('Data deleted from: %s',catrec.CaseDescription));
            end
        end
%%
        function qcData(obj,classrec,catrec,muicat) %#ok<INUSD>
            %quality control a dataset
            warndlg('No qualtiy control defined for this format');
        end    
%%
        function tabPlot(obj,src)
            %generate plot for display on Q-Plot tab
            dst = obj.Data.Dataset;
            
            if height(dst)>1
                %propmpt user to select timestep
                list = dst.DataTable.Properties.RowNames;
                irec = listdlg('PromptString','Select timestep:',...
                               'Name','Tab plot','SelectionMode','single',...
                               'ListSize',[250,100],'ListString',list);
            else
                irec = 1;
            end

            xi = dst.Dimensions.X;
            yi = dst.Dimensions.Y;
            zi = squeeze(dst.Z(irec,:,:))';
            if isempty(xi), return; end
            
            %clean up tab          
            ht = findobj(src,'Type','axes');
            delete(ht);
            %now plot results
            ax = axes('Parent',src,'Tag','Profile');
            contourf(xi,yi,zi);
            colormap('parula')
            shading flat
            caxis([-8,2]);
            colorbar
            xlabel('Length (m)'); 
            ylabel('Width (m)'); 
            title(dst.Description);
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot
        end     
        
    end
%%
    methods (Access = private)
         function setDataSetRecord(obj,muicat,dataset,datatype)
            %assign dataset to class Data property and update catalogue
            if isstruct(dataset)
                obj.Data = dataset;   %can be struct of multiple tables
            else
                obj.Data.Dataset = dataset;  
            end 
            classname = metaclass(obj).Name;            
            %add record to the catalogue and update mui.Cases.DataSets
            caserec = addRecord(muicat,classname,datatype);
            casedef = getRecord(muicat,caserec);
            obj.CaseIndex = casedef.CaseID;
            datasets = fieldnames(obj.Data);
            for i=1:length(datasets)
                if isa(obj.Data.(datasets{i}),'dstable')
                    obj.Data.(datasets{i}).Description = casedef.CaseDescription;
                end
            end
            %
            if isempty(muicat.DataSets) || ~isfield(muicat.DataSets,classname) ||...
                    isempty(muicat.DataSets.(classname))
                idrec = 1;
            else
                idrec = length(muicat.DataSets.(classname))+1;
            end
            muicat.DataSets.(classname)(idrec) = obj;           
        end  
%%
        function datasetname = getDataSetName(obj)
            %check whether there is more than one dstable and select
            dataset = 1;
            datasetnames = fieldnames(obj.Data);
            if length(datasetnames)>1
                promptxt = {'Select dataset'};
                title = 'Save dataset';
                [dataset,ok] = listdlg('PromptString',promptxt,...
                           'SelectionMode','single','Name',title,...
                           'ListSize',[300,100],'ListString',datasetnames);
                if ok<1, return; end       
            end
            datasetname = datasetnames{dataset};
        end
%%        
       function [xaxis,yaxis] = getplotaxes(~,r,x,y)
           %define axes with user defined range or a default range (N=100)
           if length(r.Xaxis)==3
               %range defined as x0:dx:xN
               xaxis = r.Xaxis(1):r.Xaxis(2):r.Xaxis(3);
           else
               %number of x-intervals defined
               nint = r.Xaxis;
               vlength = max(x)-min(x);
               xaxis = min(x):vlength/nint:max(x);
           end

           if length(r.Yaxis)==3
               %range defined as y0:dy:yN
               yaxis = r.Yaxis(1):r.Yaxis(2):r.Yaxis(3);
           else
               %number of y-intervals defined
               nint = r.Yaxis(1);
               vwidth = max(y)-min(y);
               yaxis = min(y):vwidth/nint:max(y);
           end
       end        
%%
       function griddata = formatGridData(~,data,isflip)
           %format grid to load into a dstable as single row per grid
           %saves data as z array [1,Nx,Ny] for each time step
           %X and Y assumed fixed and saved as dimensions
           X = unique(data(:,1));
           Y = unique(data(:,2));
           Nx = length(X); Ny = length(Y);
           Z = reshape(data(:,3),Nx,Ny);
           if isflip
               Z = flipud(Z);
           end
           griddata = reshape(Z,1,Nx,Ny);
           %check input plot
           %                figure;
           %                contourf(X,Y,Z');
       end
%%
        function [data,header,filename] = readInputData(~,filename,nhead)
            %load selected grid data
            fid = fopen(filename, 'r');
            if fid<0
                errordlg('Could not open file for reading','File write error','modal')
            end

            %read header and find number of columns in file
            header = cell(nhead,1);
            for i=1:nhead
                header{i} = fgetl(fid);
            end
            
            if ~contains(strtrim(header{1}),'%')
                warndlg('File does not have format specified in header');
                return;
            end
%             %should read white space, tabs and comma separated formats
%             colhead  = strsplit(header{nhead},{'\s*','\t',',\s*'},...
%                 'DelimiterType','RegularExpression','CollapseDelimiters',true);
%             ncols = length(colhead);
            ncols = 3;
            colhead = {NaN};
            %read numeric data
            if strcmpi(colhead{1},'Year')
                dataSpec = ['%{yyyy}D', repmat(' %f',1,(ncols-1))];  
            elseif strcmpi(colhead{1},'Date')
                dataSpec = ['%{dd/MM/yyyy}D', repmat(' %f',1,(ncols-1))];    
            elseif strcmpi(colhead{1},'Date/Time')
                dataSpec = ['%{dd/MM/yyyy}D %{HH:mm:ss}D', repmat(' %f',1,(ncols-1))];  
            else
                dataSpec = header{1};
            end
            data = textscan(fid,dataSpec);
            data = cell2mat(data);
            %
            fclose(fid);
        end     
        
%%
        function ax = plotGrid(~,hf,x,y,z)
            %create plot of perturpation from initial surface
            ax = axes(hf);
            pcolor(ax,x,y,z)
            shading interp
            colormap('parula')
            colorbar
            xlabel('Length (m)'); 
            ylabel('Width (m)')
        end
%%
        function [subdomain,defaults] = getSubDomain(obj,x,y,z)
            %Allow the user to interactively adjust the grid domain to sample
            subdomain0 = [min(x),max(x),min(y),max(y)];
            figtitle = sprintf('Subgrid selection');
            promptxt = 'Accept subgrid definition';
            tag = 'PlotFig'; %used for collective deletes of a group
            butnames = {'Yes','No'};
            position = [0.2,0.4,0.4,0.4];
            [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position);
            ax = plotGrid(obj,h_plt,x,y,z');
            [~,~,~,ixo,iyo] = getsubgrid(x,y,z',subdomain0);
            hold on
            plot(ax,x(ixo),y(iyo),'--r','LineWidth',0.8)
            hold off
            ok = 0; subdomain = subdomain0;
            while ok<1
                waitfor(h_but,'Tag');
                if ~ishandle(h_but) %this handles the user deleting figure window
                    ok = 1;         %continue using default subdomain
                elseif strcmp(h_but.Tag,'No')
                    %Get user to redfine subgrid
                    promptxt = {'Min X','Max X','Min Y','Max Y'};
                    title = 'Define subgrid';
                    defaults = string(subdomain');
                    asub = inputdlg(promptxt,title,1,defaults);
                    subdomain = cellfun(@str2double,asub)';
                    if any(subdomain(1,[1,3])<subdomain0(1,[1,3])) || ...
                                any(subdomain(1,[2,4])>subdomain0(1,[2,4]))
                        warndlg('Selection out of bounds. Please make a new seslection')
                        subdomain = subdomain0;
                    end
                    [~,~,~,ixo,iyo] = getsubgrid(x,y,z',subdomain);
                    h_sd = findobj(ax,'Type','line');  %remove existing subdomain rectangle
                    delete(h_sd);
                    hold on
                    plot(ax,x(ixo),y(iyo),'--r','LineWidth',0.8) %updated subdomain bounding rectangle
                    hold off
                    h_but.Tag = '';
                else
                    ok = 1;
                end
            end
        end
%%        
        function dsp = setDSproperties(~)
            %define the metadata properties for the demo data set
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...
                'Name',{'Z'},...                  
                'Description',{'Elevation'},...
                'Unit',{'mAD'},...
                'Label',{'Elevation (mAD)'},...
                'QCflag',{'raw'});  
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'y'},...
                'Label',{'Time (yr)'},...
                'Format',{'y'});       
            dsp.Dimensions = struct(...    
                'Name',{'X','Y'},...
                'Description',{'X-axis','Y-axis'},...
                'Unit',{'m','m'},...
                'Label',{'X-axis (m)','Y-axis (m)'},...
                'Format',{'',''}); 
        end
    end
end