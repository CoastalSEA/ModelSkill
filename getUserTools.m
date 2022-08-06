function getUserTools(mobj)
%
%-------function help------------------------------------------------------
% NAME
%   getUserTools.m
% PURPOSE
%   user functions to do additional analysis on data loaded in ModelSkill
% USAGE
%   getUserTools(mobj)
% INPUTS
%   mobj - ModelUI instance
% OUTPUT
%   Network Analysis: 
%   Rhythmic Forms: 
% NOTES
%   uses tools in GDinterface 
% SEE ALSO
%   getInletTools. Called in ModelSkill
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2022
%--------------------------------------------------------------------------
%
    ok = 1;
    while ok>0
        usertools = {'Network Analysis','Rhythmic Forms'};                
        [selection,ok] = listdlg('Liststring',usertools,...
                                 'PromptString','Select tool:',...  
                                 'ListSize',[180,100],...
                                 'SelectionMode','single');
        if ok==0, continue; end  %user cancelled selection
        src = usertools{selection};
        %
        switch src
            case 'Network Analysis'
                getNetworkStats(mobj); 
            case 'Rhythmic Forms'
                getRhythmicForm(mobj);
        end
    end
end
%%
function getNetworkStats(mobj)
    %select case and call network_count function
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to analyse:','Select timestep:'};
    [cobj,~,irow] = selectCaseDatasetRow(mobj.Cases,[],gridclasses,promptxt,1);
    if isempty(cobj) || isempty(irow), return; end

    casedesc = cobj.Data.Form.Description;
    timetxt = cobj.Data.Form.DataTable.Properties.RowNames{irow};
    casedesc = sprintf('%s at %s',casedesc,timetxt);
    [options,casevars] = getNetworkVars(casedesc);
    bathy = squeeze(cobj.Data.Form.Z(irow,:,:));
    network_count(bathy,options,casevars);
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
    casevars.X_centre = '85';
    casevars.Y_centre = '1';
    casevars.min_rad = '12';
    casevars.rad_int = '2';
    casevars.grid_size = '100';
    casevars.HW = '1';
    casevars.isplot = '0';
    casevars.desc = casedesc;

    defaultvals = struct2cell(casevars);
    promptxt = {'X-index for circle centre','Y-index for circle centre',...
                'Minimum radius','Radius interval',...
                'Grid size (m)','High water level (mOD)',...
                'Include all plots (0/1)','Case description'};
    answer = inputdlg(promptxt,'Case variables',1,defaultvals);   
    cvar = cellfun(@str2double,answer(1:end-1),'UniformOutput',false);
    cvar{end} = logical(cvar{end});
    cvar = [cvar;answer(end)];
    fields = fieldnames(casevars);
    casevars = cell2struct(cvar,fields);
end
%%
function hf = getRhythmicForm(mobj)
    %plot rhythmic forms on a horizontal surface (ie bed perturbation)
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to analyse:','Select timestep:'};
    [cobj,~,irow] = selectCaseDatasetRow(mobj.Cases,[],gridclasses,promptxt,1);
    if isempty(cobj) || isempty(irow), return; end
    grid = getGrid(cobj,irow);
    
    casedesc = cobj.Data.Form.Description;
    timetxt = cobj.Data.Form.DataTable.Properties.RowNames{irow};
    casedesc = sprintf('%s at %s',casedesc,timetxt);

    %find the highest and lowest elevation in the y coordinate to define
    %initial slope and return differences from initial slope
    [X,Y,zdiff] = getGridPerturbation(grid);
    %retrieve description from Results Case
    %now examine rhythmic spacing
    SD = getRhythmicSpacing(cobj,grid,zdiff,casedesc);
    if isempty(SD), return; end
    %get the volume change relative to the zero horizontal surface
    [posvol,negvol] = getVolumeChange(grid,zdiff);
    %create plot of perturbation from initial surface     
    hf = rhythmicFormPlot(X,Y,zdiff);
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
function SD = getRhythmicSpacing(obj,grid,zdiff,casedesc)
    %compute the rhythmic form spacing and amplitude
    SD = [];
    %create subgrid selection plot
    [subdomain,sublimitxt] = getSubDomain(obj,grid);
    if isempty(subdomain), return; end
    casedesc = sprintf('Subgrid of %s\nUsing %s',casedesc,sublimitxt);

    grid.z = zdiff';
    [subgrid,SD.x,SD.y] = getsubgrid(grid,subdomain);
    subz = subgrid.z';
    
    %compute stats for each row over the ysub domain
    pamp = zeros(size(subgrid.y)); wavelength = pamp; stdpamp = pamp;
    namp = pamp; stdnamp = pamp;    
    nrow = length(subgrid.y);
    for i=1:nrow
        var = subz(i,:);
        [pks,~] = getRowPeaks(var);
        %compute the mean range and wavelength for each row
        pamp(i) = mean(pks(pks>0),'omitnan');
        stdpamp(i)  = std(pks(pks>0),'omitnan');
        namp(i) = mean(pks(pks<0),'omitnan');
        stdnamp(i)  = std(pks(pks<0),'omitnan');
        wavelength(i) = (subgrid.x(end)-subgrid.x(1))/length(pks)*2;
    end

    %amplitude and wavelength plot
    hf = figure('Name','Rhythmic properties', ...
                'Units','normalized', ...
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    yyaxis left
    plot(ax,subgrid.y,pamp,'-','Color',mcolor(3),'DisplayName','Positive amplitude') %yellow
    hold on 
    plot(ax,subgrid.y,pamp+stdpamp,':','Color',mcolor(3),'DisplayName','+/- 1st.dev.') %yellow
    p1 = plot(ax,subgrid.y,pamp-stdpamp,':','Color',mcolor(3)); %yellow
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    plot(ax,subgrid.y,namp,'-','Color',mcolor(1),'DisplayName','Negative amplitude') %dark blue
    plot(ax,subgrid.y,namp+stdnamp,':','Color',mcolor(1),'DisplayName','+/- 1st.dev.')%dark blue
    p1 = plot(ax,subgrid.y,namp-stdnamp,':','Color',mcolor(1));%dark blue
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    hold off
    xlabel('Width (m)');
    ylabel('Mean amplitude and +/-St.Dev. (m)');
    yyaxis right
    plot(ax,subgrid.y,wavelength,'.-.','DisplayName','Wavelength')
    ylabel('Mean wavelength (m)');  

    m1 = sprintf('Subdomain wavelength and amplitude for %s',casedesc);
    m2 = sprintf('Wavelength: mean=%.1f; mode=%0.1f; min=%0.1f.',...
                    mean(wavelength), mode(wavelength),min(wavelength));
    m3 = sprintf('Positive mean amplitude = %.2f; Negative mean amplitude = %.2f',...
                               mean(pamp,'omitnan'),mean(namp,'omitnan'));
    legend('Location','southwest')
    titletxt = sprintf('%s\n%s\n%s',m1,m2,m3);
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
function hf = rhythmicFormPlot(X,Y,Z)
    %create plot of perturpation from initial surface
    hf = figure('Name','Rhythmic properties', ...
                'Units','normalized', ...
                'Resize','on','HandleVisibility','on', ...
                'Tag','PlotFig');
    ax = axes(hf);
    pcolor(ax,X,Y,Z)
    shading interp
    colormap('parula')
    colorbar
    % axis('equal')
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
    pks = zeros(size(first));
    for j= 1:length(first)
        if rlen(j)>1
            firstpk = var(locs(first(j),1));
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


