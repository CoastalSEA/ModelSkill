function getUserTools(mobj,src,caseid)
    %user functions to do additional analysis on data loaded in ModelSkill
    %mobj - ModelUI instance
    %src - usertools selection defined in call
    %useCase - selected grid to add properties to
    if nargin<2 || isempty(src)
        usertools = {'Gross properties','Tabulate properties','Plot properties'...
                     'Plot hypsometry','Network analysis','Cusp plot'};
                     
        [selection,ok] = listdlg('Liststring',usertools,...
                                 'PromptString','Select tool:',...  
                                 'ListSize',[180,100],...
                                 'SelectionMode','single');
        if ok==0, return; end  %user cancelled selection
        src = usertools{selection};
        caseid = [];
    end
    switch src
        case 'Gross properties'
            getGrossProperties(mobj,caseid);
        case 'Tabulate properties'
            getPropsTable(mobj,true);
        case 'Plot properties'
            getPropsPlot(mobj);
        case 'Plot hypsometry'
            getHypsPlot(mobj,caseid);
        case 'Network analysis'
            getNetworkStats(mobj,caseid); 
        case 'Cusp plot'
            getCuspPlot(mobj,caseid);
    end
end
%%
function getGrossProperties(mobj,caseid)
    %prompt user to select an existing form model to modify
    muicat = mobj.Cases;
    if isempty(caseid)
        promptxt = 'Select a Model:';
        [caserec,~] = selectCase(muicat,promptxt,'single',0);     
        if isempty(caserec), return; end
    else
        caserec = caseRec(mobj.Cases,caseid);
    end
    wl = getTidalLevels();    
    if isempty(wl), return; end
    wls = struct2cell(wl)';
    
    [cobj,classrec,~] = getCase(muicat,caserec);
    dst = cobj.Data.Dataset;
    nrec = height(dst);
    grid.x = dst.Dimensions.X;
    grid.y = dst.Dimensions.Y;
    allz = dst.Z;
%     hyps{1,4} = 0;
%     gprop{1,15} = 0;
    wldst = dstable(wls{:},'DSproperties',PropsDSP('Levels'));
    grid.z = squeeze(allz(1,:,:));
    hyps = getHypsometry(grid,wl);
    hypsdst = dstable(hyps{:},'DSproperties',PropsDSP('Hyps')); 
    gprop = struct2cell(grossProperties(grid,wl,hyps))';    
    gpropdst = dstable(gprop{:},'DSproperties',PropsDSP('Gross'));
    
    
    for i=2:nrec  %process all rows in the dstables
        wldst.DataTable = [wldst.DataTable;wls];
        grid.z = squeeze(allz(i,:,:));
        hyps = getHypsometry(grid,wl);
        hypsdst.DataTable = [hypsdst.DataTable;hyps];
        gprop = struct2cell(grossProperties(grid,wl,hyps))';
        gpropdst.DataTable = [gpropdst.DataTable;gprop];
    end
    
    %add RowNames to each dstable
    wldst.RowNames = dst.RowNames;
    hypsdst.RowNames = dst.RowNames;
    gpropdst.RowNames = dst.RowNames;
  
    cobj.Data.WaterLevels = wldst;          %assign hypsometry results
    cobj.Data.GrossProps = gpropdst;        %asign gross property results
    cobj.Data.Hypsometry = hypsdst;         %assign hypsometry results

    updateCase(muicat,cobj,classrec,'Gross properties added');
end
%%
function getNetworkStats(mobj,caseid)
    if isempty(caseid)
        caserec = [];
    else
        caserec = find(caseid==mobj.Cases.CaseID);      
    end
    [grid,timetxt,caserec,~] = MS_GridData.getCaseGridData(mobj.Cases,caserec);
    if isempty(grid), return; end

    casedesc = mobj.Cases.Catalogue.CaseDescription(caserec);
    casedesc = sprintf('%s at %s',casedesc,timetxt);
    [options,casevars] = getNetworkVars(casedesc);
    network_count(grid.z,options,casevars);
end
%%
function [options,casevars] = getNetworkVars(casedesc)
    %define analysis parameters
    % options - struct to define analysis parameters used in geonet_diffusion
    %   method -'lin':  Linear diffusion (constant c=1).
    %           'pm1': perona-malik, c=exp{-(|grad(J)|/K)^2} [PM90]
    %           'pm2': perona-malik, c=1/{1+(|grad(J)|/K)^2} [PM90]
    %           'tukey: ??
    %           'rmp': complex valued - ramp preserving [GSZ01]
    %   lambda - edge threshold parameter
    %   nint - number of iterations
    %   dt - time increment (0 < dt <= 0.25, default 0.2)
    %   sigma2 - if present, calculates gradients of diffusion coefficient
    %          convolved with gaussian of var sigma2 (Catte et al [CLMC92])J
    %
    % casevars - struct to define the parameters used for channel counting
    %   X_centre -grid index to centre of circle along x-axis
    %   Y_centre -grid index to centre of circle along y-axis
    %   min_rad - minimum radius (in grid intervals)
    %   rad_int - number of grid intervals used to generate radii
    %   grid_size - grid dimension (m)
    %   HW - high water level (m)
    %   casedesc - description of selected case

    % assign values to options structure for use in geonet_diffusion
    options.method = 'pm2';
    options.lambda = '1';
    options.nint = '2';
    options.dt = '0.2';
    options.sigma = '0';

    defaultvals = struct2cell(options);
    promptxt = {'Method','Lambda','No of intervals','Time increment','Sigma'};       
    answer = inputdlg(promptxt,'Network extraction',1,defaultvals);
    if isempty(answer), return; end

    %assign variables based on user definition
    opts = answer(2:end);
    opts = cellfun(@str2double,opts,'UniformOutput',false);
    opts = [answer(1);opts];
    fields = fieldnames(options);
    options = cell2struct(opts,fields);

    % assign values to casevars structure for use in channel counting
    casevars.X_centre = '0';
    casevars.Y_centre = '0';
    casevars.min_rad = '12';
    casevars.rad_int = '2';
    casevars.grid_size = '100';
    casevars.HW = '1';
    casevars.desc = casedesc;

    defaultvals = struct2cell(casevars);
    promptxt = {'X-index for circle centre','Y-index for circle centre',...
                'Minimum radius','Radius interval',...
                'Grid size (m)','High water level (mOD)','Case description'};
    answer = inputdlg(promptxt,'Case variables',1,defaultvals);   
    cvar = cellfun(@str2double,answer(1:end-1),'UniformOutput',false);
    cvar = [cvar;answer(end)];
    fields = fieldnames(casevars);
    casevars = cell2struct(cvar,fields);
end
%%
function wl = getTidalLevels()
    %prompt user for levels to used to compute volumes and areas at high
    %and low water
    prmptxt = {'High water level:','Mean tide level:','Low water level:'};
    dlgtitle = 'Water levels';
    defaultvals = {'1','0','-1'};
    answer = inputdlg(prmptxt,dlgtitle,1,defaultvals);
    if isempty(answer)
        wl = [];  
        return;
    end
    wl.HW = str2double(answer{1});
    wl.MT = str2double(answer{2});
    wl.LW = str2double(answer{3});
end
%%
function hyps = getHypsometry(grid,wl)
    %compute area hypsometry from grid data
    histint = 0.1;
    lowlim = -10.0;
    delx = grid.x(2)-grid.x(1);
    dely = grid.y(2)-grid.y(1);
    % this is from ChannelForm model where wl varies along channel
    %     % range for histogram data - general
    %     lowlim = floor(min(min(grid.z)));
    %     %the upper limit is problematical when hw level varies
    %     %using the maximum means some area above hw is included
    %     uplim = ceil(max(wl.HW));        
    %     zedge = lowlim:histint:uplim;    
    %     zedge(zedge>(wl.HW+histint)) = [];
    zedge = lowlim:histint:wl.HW;
    % calculate histogram and format output
    zed  = reshape(grid.z,numel(grid.z),1);
    % zed  = zed(zed<uplim);     %control range for hypsometry histogram
    % zed  = zed(zed>=lowlim);
    zhist = histcounts(zed,zedge); %bin counts for each interval defined by zedge
    zhist = zhist*delx*dely;     %scale occurences by grid cell area
    zsurf = cumsum(zhist);
    zvol = zsurf*histint;        %volume of each slice
    zvol = cumsum(zvol);         %cumulative volume below each elevation
    zcentre = movsum(zedge,[0,1])/2;
    zcentre(end) = [];
    hyps  = {zcentre,zhist,zsurf,zvol};    
end
%%
function grossprops = grossProperties(grid,wl,hyps)
    %compute the gross properties of the selected bathymetry
    histint = 0.1+eps;
    am = (wl.HW-wl.LW)/2;         %tidal amplitude(m)
    z0 = wl.MT;         %mean tide level(m)
    %lowlim = inp.zm; uplim = inp.zhw;%define range to be entire model domain
    zcentre = hyps{:,1};
    zsurf = hyps{:,3};
    zvol = hyps{:,4};

    idh = find(zcentre<(wl.HW+histint/2) & zcentre>=(wl.HW-histint/2),1,'first');
    idl = find(zcentre<=(wl.LW+histint/2) & zcentre>=(wl.LW-histint/2),1,'first');
    ido = find(zcentre<=(z0+histint/2) & zcentre>=(z0-histint/2),1,'first');

    r.Shw = zsurf(idh);               %surface area at high water
    r.Slw = zsurf(idl);               %surface area at low water
    r.Vhw = zvol(idh);                %volume at high water
    r.Vlw = zvol(idl);                %volume at low water
    r.Pr  = r.Vhw-r.Vlw;              %volume of tidal prism
    beta = 1;                         %assume Schw/Sclw = 1
    r.Gam = (r.Slw/r.Shw)^3*(r.Vhw/r.Vlw)^2*beta;  %Dronkers gamma (~1)
    r.Vs  = (r.Vhw-r.Vlw)-2*am*r.Slw; %storage volume over intertidal
    r.Vc  = r.Vlw+2*am*r.Slw;         %channel volume
    r.VsoVc = r.Vs/r.Vc;              %ratio fo storage to channel volumes
    r.a   = am;                       %tidal amplitude
    r.h = zvol(ido)/zsurf(ido);       %hydraulic depth at mtl
    r.aoh = am/r.h;                   %tidal amplitude/hydraulic depth    
    x0 = min(grid.x);
    [r.W0,r.A0] = getCSAat_z0_X(grid,z0,x0); %width and csa at mouth
    r.PoA = r.Pr/r.A0;                %Prism to CSA ratio
    grossprops = r;                   %assign structure to output
end
%%
function [w,csa] = getCSAat_z0_X(grid,z0,x)
    %get thw width and cross-sectional area at distance x and elevation z0
    yi = grid.y;
    dely = yi(2)-yi(1);  %grid interval
    idx = (x==grid.x);
    zi = grid.z(idx,:);
    zi(zi>z0) = NaN;
    w = sum(~isnan(zi)*dely);
    dwl = z0 - zi;
    csa = trapz(yi(~isnan(dwl)),dwl(~isnan(dwl)));
end
%%
function getHypsPlot(mobj,caseid)
    %plot the grid hypsometry
    if isempty(caseid)
        caserec = [];
    else
        caserec = find(caseid==mobj.Cases.CaseID);      
    end
    [grid,timetxt,caserec,irec] = MS_GridData.getCaseGridData(mobj.Cases,caserec);
    if isempty(grid), return; end
    %retrieve additional data from data set
    cobj = getCase(mobj.Cases,caserec);
    wl = cobj.Data.WaterLevels.DataTable(irec,:);
    hyps = cobj.Data.Hypsometry.DataTable(irec,:);    
    zcentre = hyps.zc;
    zcentre(end,1) = zcentre(end,1)+0.1; %crude offset to highest point so that HW is visible
    zsurf = hyps.zs;
    zvol = hyps.zv; 
    %retrieve description from Results Case
    casedesc = mobj.Cases.Catalogue.CaseDescription(caserec);
    casedesc = sprintf('%s at %s',casedesc,timetxt);
    
    figure;
    plot(zvol,zcentre,'-r');
    hold on
    plot(zsurf,zcentre,'-.b');
    
    %add water levels at mouth
    plot(xlim, wl.HW*[1 1],':','Color',[0.7,0.7,0.7]);
    plot(xlim, wl.MT*[1 1],'--','Color',[0.8,0.8,0.8]);
    plot(xlim, wl.LW*[1 1],':','Color',[0.7,0.7,0.7]);
    
    xlabel('Volume (m^3) and Area (m^2)'); 
    ylabel('Elevation (mAD)');
    L1 = sprintf('Vol - %s',casedesc);
    L2 = sprintf('Area - %s',casedesc);
    legend(L1,L2,'Location','SouthEast');
    title('System hypsometry','FontWeight','normal','FontSize',10);
    hold off
end
%%
function tabledata = getPropsTable(mobj,istab)
    %tabulate the gross properties in a figure
    rN = 15; %No. of rows in the gross properties struct
    muicat = mobj.Cases;
    promptxt = 'Select a Model:';
    [caserec,~] = selectCase(muicat,promptxt,'single',0);     

    %retrieve additional data from data set
    cobj = getCase(mobj.Cases,caserec);
    if isfield(cobj.Data,'GrossProps')
        tabledata = cobj.Data.GrossProps;
    else
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
        tabledata = getPropsTable(mobj,false); 
        if isempty(tabledata)
            hold off
            legend
            return; 
        end
        PromptText = 'Select Variable'; %default if no varargin
        ListSize = [300,250];           %default if no varargin
        [varnum,ok] = listdlg('Name','Variable', ...
                'ListSize',ListSize,...
                'PromptString',PromptText, ...
                'SelectionMode','single', ...
                'ListString',tabledata.VariableNames);  
        if ok~=0
            data = tabledata.DataTable{:,varnum};
            x = tabledata.RowNames;
            xlabel('Time step (years)')
            legtext = inputdlg('Legend text','Plot properties',1,{num2str(count)});
            plot(ax,x,data,'LineWidth',1,'DisplayName',legtext{1});
            ylabel(tabledata.VariableLabels{varnum});
            answer = questdlg('Add another dataset?','Plot properties','Yes');
            if strcmp(answer,'No')
                select = 0;
            else
                count = count+1;
            end
        end
    end
    hold off
    legend
end
%%
function hf = getCuspPlot(mobj,caseid)
    %plot a cusp model on a horizontal surface (ie bed perturbation)
    if isempty(caseid)
        caserec = [];
    else
        caserec = find(caseid==mobj.Cases.CaseID);      
    end
    [grid,timetxt,caserec,~] = MS_GridData.getCaseGridData(mobj.Cases,caserec);
    if isempty(grid), return; end
    %retrieve description from Results Case
    casedesc = mobj.Cases.Catalogue.CaseDescription(caserec);
    casedesc = sprintf('%s at %s',casedesc,timetxt);
    %find the highest and lowest elevation in the y coordinate to define
    %initial slope and return differences from initial slope
    [X,Y,zdiff] = getGridPerturbation(grid);
    %retrieve description from Results Case
    %now examin cusp spacing
    SD = getCuspSpacing(mobj,grid,zdiff,casedesc);
    if isempty(SD), return; end
    %get the volume change relative to the zero horizontal surface
    [posvol,negvol] = getVolumeChange(grid,zdiff);
    %create plot of perturbation from initial surface     
    hf = cuspPlot(X,Y,zdiff);
    ax = findobj(hf,'Type','axes');
    hold on
    plot(ax,grid.x(SD.x),grid.y(SD.y),'--r')
    hold off    
    m1 = sprintf('Perturbation from initial plane for %s',casedesc);
    m2 = sprintf('Positive volume change = %.1f; Negative volume change = %.1f',...
                            posvol,negvol);
    m3 = 'Dashed lines: subdomain used to estimate amplitude and wavelength';                    
    titletxt = sprintf('%s\n%s\n%s',m1,m2,m3);
    title(titletxt)    
end
%%
function SD = getCuspSpacing(mobj,grid,zdiff,casedesc)
    %compute the cusp spacing and amplitude
    %extract Skill parameters from MS_RunParameters class
    SD = [];
    msgtxt = 'Averaging window has not been defined';
    if isfield(mobj.Inputs,'MS_RunParams')
        subdomain = mobj.Inputs.MS_RunParams.skillsubdomain;
        if length(subdomain)~=4
            warndlg(msgtxt); return;
        end        
    else
        warndlg(msgtxt);  return
    end  
    [xsub,ysub,zsub,SD.x,SD.y] = getsubgrid(grid.x,grid.y,zdiff,subdomain);
    
    %compute stats for each row over the ysub domain
    pamp = zeros(size(ysub)); wavelength = pamp; stdpamp = pamp;
    namp = pamp; stdnamp = pamp;    
    nrow = length(ysub);
    for i=1:nrow
        var = zsub(i,:);
        [pks,~] = getRowPeaks(var);
        %compute the mean range and wavelength for each row
        pamp(i) = mean(pks(pks>0),'omitnan');
        stdpamp(i)  = std(pks(pks>0),'omitnan');
        namp(i) = mean(pks(pks<0),'omitnan');
        stdnamp(i)  = std(pks(pks<0),'omitnan');
        wavelength(i) = (xsub(end)-xsub(1))/length(pks)*2;
    end

    %amplitude and wavelength plot
    hf = figure('Name','Cusp properties', ...
                'Units','normalized', ...
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    yyaxis left
    plot(ax,ysub,pamp,'-','Color',[0.929 0.694 0.125])
    hold on 
    plot(ax,ysub,pamp+stdpamp,':','Color',[0.929 0.694 0.125])
     plot(ax,ysub,pamp-stdpamp,':','Color',[0.929 0.694 0.125])
    plot(ax,ysub,namp,'-','Color',[0 0.447 0.741])
    plot(ax,ysub,namp+stdnamp,':','Color',[0 0.447 0.741])
    plot(ax,ysub,namp-stdnamp,':','Color',[0 0.447 0.741])
    hold off
    xlabel('Width (m)');
    ylabel('Mean amplitude and +/-St.Dev. (m)');
    yyaxis right
    plot(ax,ysub,wavelength,'.-.')
    ylabel('Mean wavelength (m)');  
%     casedesc = mobj.Cases.CaseDescription{useCase};
    m1 = sprintf('Subdomain wavelength and amplitude for %s',casedesc);
    m2 = sprintf('Wavelength: mean=%.1f; mode=%0.1f; min=%0.1f.',...
                    mean(wavelength), mode(wavelength),min(wavelength));
    m3 = sprintf('Positive amplitude = %.2f; Negative amplitude = %.2f',mean(pamp),mean(namp));
    m4 = 'Lines are +ve/-ve amplitude and dash-dot is wavelength';
    titletxt = sprintf('%s\n%s\n%s\n%s',m1,m2,m3,m4);
    title(titletxt)
end
%%
function [X,Y,zdiff] = getGridPerturbation(grid)
    %compute the variabiton of the surface from the initial plane
    %find the highest and lowest elevation in the y coordinate
    zmin = min(min(grid.z));
    zmax = max(max(grid.z));
    nlen = length(grid.x);
    y_offset = min(grid.y);
    ylength = max(grid.y)-y_offset;
    
    zfun = @(y) (zmax-zmin)*(y-y_offset)/ylength;
    zflat = repmat(zfun(grid.y),1,nlen);
    zdiff = grid.z'-zflat;
    [X,Y] = meshgrid(grid.x,grid.y);
end
%%
function [posvol,negvol] = getVolumeChange(grid,zdiff)
    %compute the change in volume relative to the zero surface    
    zpos   = zdiff.*(zdiff>0);
    posvy  = trapz(grid.y,zpos,1);
    posvol = trapz(grid.x,posvy);
    zneg   = zdiff.*(zdiff<0);
    negvy  = trapz(grid.y,zneg,1);
    negvol = trapz(grid.x,negvy);
end
%%
function hf = cuspPlot(X,Y,Z)
    %create plot of perturpation from initial surface
    hf = figure('Name','Cusp properties', ...
                'Units','normalized', ...
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    pcolor(ax,X,Y,Z)
    shading interp
    colormap('parula')
    colorbar
%     axis('equal')
    xlabel('Length (m)'); 
    ylabel('Width (m)')
end
%%
function [pks,idx] = getRowPeaks(var)
    %find alternate peaks and troughs and return peaks and index of peak
    [plocs,~] = peakseek(var,1,0.01); %(variable, min distance, min amplitude)
    [nlocs,~] = peakseek(-var,1,0.01);
    %create array of peak locations and whether they are +ve or -ve
    locs = sortrows([[plocs,nlocs]',[true(size(plocs)),false(size(nlocs))]']);
    %find run lengths (ie first occurrence of several +ve/-ve peaks)
    %(Code from https://uk.mathworks.com/matlabcentral/answers/469659-count-the-adjacent-same-elements-in-a-vector)
    A = locs(:,2);
    first = find([true diff(A')~=0]);
    rlen = diff([first numel(A)+1]); % run-lengths, elements given by C = A(first); 
    %find the peaks for each alternate max/min in the row
    pks = zeros(size(first)); xpos = pks;
    for j= 1:length(first)
        if rlen(j)>1
            firstpk = var(locs(first(j),1));
%                 idfirst = locs(first(j),1);
            idlast = first(j)+rlen(j)-1;
            if firstpk>=0
                [~,idm] = max(var(locs(first(j):idlast)));
            else
                [~,idm] = min(var(locs(first(j):idlast)));
            end
            idx = locs(first(j)+idm-1);
        else               
            idx = locs(first(j),1);
        end
        pks(j) = var(idx);
        % xpos(j) = xsub(idx);        
%         if j>1
%             range(j) = abs(pks(j)-pks(j-1));
%         else
%             range(j) = NaN;
%         end
    end
    %check plot
    % figure
    % plot(xsub,var)
    % hold on
    % plot(xsub(plocs),ppks,'.b')
    % plot(xsub(nlocs),-npks,'xr')
    % plot(xpos,pks,'og')
    % hold off
end
%%
function dsp = PropsDSP(type)
            %define the metadata properties for the demo data set
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]);  
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'y'},...
                'Label',{'Time (yr)'},...
                'Format',{'y'});       
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{')'},...
                'Format',{''}); 
            switch type
                case 'Gross'
                    dsp.Variables = struct(...
                        'Name',{'Shw','Slw','Vhw','Vlw','Pr','Gam',...
                                'Vs','Vc','VsoVc','a','h','aoh','W0','CSA0','PoA'},...                  
                        'Description',{'Surface area at high water',...
                                'Surface area at low water',...
                                'Volume at high water',...
                                'Volume at low water',...
                                'Volume of tidal prism',...
                                'Dronkers gamma',...
                                'Storage volume over intertidal',...
                                'Channel volume to mtl',...
                                'Ratio of storage volume to channel volume',...
                                'Tidal amplitude',...
                                'Hydraulic depth at mtl',...
                                'Tidal amplitude to hydraulic depth ratio',...
                                'Width at mouth',...
                                'CSA at mouth',...
                                'Prism to CSA ratio'},...
                        'Unit',{'m^2','m^2','m^3','m^3','m^3','-','m^3',...
                                'm^3','-','m','m','-','m','m^2','m'},...
                        'Label',{'Surface area (m^2)','Surface area at (m^2)',...
                                'Volume (m^3)','Volume (m^3)','Volume (m^3)',...
                                'Dronkers gamma (-)','Volume (m^3)','Volume (m^3)',...
                                'Ratio of storage volume to channel volume',...
                                'Tidal amplitude (m)','Hydraulic depth (m)',...
                                'Tidal amplitude to hydraulic depth ratio',...
                                'Width (m)','Cross-sectiona area (m^2)',...
                                'Prism to CSA ratio'},...
                        'QCflag',repmat({'model'},1,15));    
                case 'Hyps'
                    dsp.Variables = struct(...
                        'Name',{'zc','zh','zs','zv'},...                  
                        'Description',{'Elevation','Area Occurrence',...
                        'Cumulative area','Cumulative volume'},...
                        'Unit',{'mAD','-','m^2','m^3'},...
                        'Label',{'Elevation (mAD)','Occurrence probability',...
                        'Area (m^2)','Volume (m^3)'},...
                        'QCflag',repmat({'model'},1,4));  
                case 'Levels'
                    dsp.Variables = struct(...
                        'Name',{'HW','MT','LW'},...                  
                        'Description',{'High water','Mean tide level','Low water'},...
                        'Unit',repmat({'mAD'},1,3),...
                        'Label',repmat({'Elevation (mAD)'},1,3),...
                        'QCflag',{'data'}); 
            end            
end



