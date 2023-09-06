function outable = network_count(bathy,options,casevars)
%
%-------function help------------------------------------------------------
% NAME
%   network_count.m
% PURPOSE
%   Extract channels from a bathymetry and perform some network analysis
%   function derived from GeoNet_phd_coundtin_channels_Ian.m provided by
%   Barend van Maanen (Apr 2020)
% USAGE
%   outable = network_count(bathy,options,casevars)
% INPUTS
%   bathy - elevations for cartesian grid
%   options - struct to define analysis parameters used in geonet_diffusion (see below)
%   casevars - struct to define the parameters used for channel counting (see below)
% OUTPUT
%   outable - table with fields Case, PeakCount, PeakDistance, AvDrainDist,
%             StdDrainDist, SlopeDrainDist, r2DrainDist
% NOTES
%   requires Matlab Image Processing Toolbox and Statistics and Machine Learning Toolbox
% SEE ALSO
%   getUserTools. Called in ModelSkill
%
% Author: Ian Townend
% CoastalSEA (c) Mar 2022, output added Aug 2023 
%--------------------------------------------------------------------------
%
%Reference: Sangireddy et al, 2016, GeoNet: An open source software for
%the automatic and objective extraction of channel heads, channel 
%network, and channel morphology from high resolution topography data, 
%Environmental Modeling and Software, 83, 58-73, doi:10.1016/j.envsoft.2016.04.026.
%
% bathy_0 - an m x n array of elevations of the bathymetry to be analysed
% tr - the tidal range (m)
% options - struct to define analysis parameters
%   method -'lin':  Linear diffusion (constant c=1).
%           'pm1': perona-malik, c=exp{-(|grad(J)|/K)^2} [PM90]
%           'pm2': perona-malik, c=1/{1+(|grad(J)|/K)^2} [PM90]
%           'tukey: ??
%           'rmp': complex valued - ramp preserving [GSZ01]
%   lambda - edge threshold parameter
%   nint - number of iterations
%   dt - time increment (0 < dt <= 0.25, default 0.2)
%   del - grid spacing (m)
%   sigma2 - if present, calculates gradients of diffusion coefficient
%          convolved with gaussian of var sigma2 (Catte et al [CLMC92])J
% casevars - struct to define the parameters used for channel counting
%   X_centre -grid index to centre of circle along x-axis
%   Y_centre -grid index to centre of circle along y-axis
%   min_rad - minimum radius (in grid intervals)
%   rad_int - number of grid intervals used to generate radii
%   grid_size - grid dimension (m)
%   HW - high water level (m)
%   isplot - generates all ouput plots if true
%   casedesc - description of selected cas 
%
%Perform Perona-Malik nonlinear filtering on the original data set 
%to obtain the regularized data set
%get the diffusion in the grid using GeoNet function diffusiion
%  
%geonet_diffusion uses imgaussfilt from Matlab Image Processing Toolbox
%geonet_qqplot uses norminv from Statistics and Machine Learning Toolbox

    v = options;
    K = geonet_diffusion(bathy, v.method, v.nint, v.lambda,...
                         v.dt, v.sigma);

                     
                     %NOT SURE WHAT THIS WAS SUPPOSED TO DO
%     [x_long, y_long]=size(bathy_0);
%     if (y_long > x_long)
%         z = -9999*ones(y_long - 1 - (x_long - 1), y_long - 1 + 1);
%         K = [K; z];
%     else
%         z = -9999*ones(x_long - 1 + 1, x_long - 1 - (y_long - 1));
%         K = [K, z];
%     end

    %Computation of the curvature on the regularized data set
    Co = curv(K,casevars.grid_size); 
    % if casevars.isplot  %check plot 
    %     figure('Name','Network','Units','normalized','Tag','PlotFig');
    %     pcolor(Co)
    %     shading flat; 
    %     colorbar; 
    %     view([0 -90])
    %     title('Curvature of the regularized data, Co')
    % end

    %Computation of statistics of curvature
    [x_Co,y_Co] = size(Co); 
    Cs_bis = reshape(Co,1,[]);
    index = Cs_bis == 0;
    Cs_bis(index) = NaN;
    jj=isnan(bathy); %bathymetry uses NaN values to represent land
    Cs_bis(jj) = NaN;
    Cs_bis = reshape(Cs_bis,x_Co,y_Co);
    
    if casevars.isplot
        figure('Name','Network','Units','normalized','Tag','PlotFig');
        pcolor(Cs_bis)
        shading flat; 
        colorbar; 
        view([0 -90])
        max_Co=max(max(Co));
        min_Co=min(min(Co));
        set(gca,'clim',[min_Co,max_Co])
        title('Curvature of the regularized data, Cs')
    end
 
    above_high_tide=bathy>casevars.HW; %grid points above the high water level
    bathy_channels = get_channels(Cs_bis,above_high_tide,casevars);
    
    results = get_number_channels(bathy_channels,casevars);
    if isempty(results)
        warndlg(sprintf('No results to plot.\nTry changing parameters used for channel counting'));
        return; 
    end
    outable = plot_network_figures(bathy,bathy_channels,results,casevars);
end
%%
function channels = get_channels(Cs_bis,above_high_tide,casevars)
    %Extract channels
    %
    %Quantile-quantile plot of curvature to set curvature threshold
    [xx,yy,qx,qy,mx,my] = geonet_qqplot(Cs_bis(~isnan(Cs_bis))); %note: nargin==1
    
    %Plotting quantile-quantile plot
    h_fig = figure;    
    plot(xx,yy,':',qx,qy,'-o',mx,my,'-.');  
    
    xlabel('Standard Normal Quantiles')
    ylabel('Quantiles of Input Sample')
    title ('QQ Plot of Sample Data versus Standard Normal')    
    
    %define a threshold (potentially based on the quantile-quantile plot)
    promptxt = sprintf('Define threshold to use\nDefault is 75 percentile of x (right circle)');
    answer = inputdlg(promptxt,'Channel extraction',1,{num2str(qx(2))});
    if isempty(answer)        %use default value if user cancels        
        threshold_1 = qx(2);  %75 percentile value
    else                      %apply user defined value
        threshold_1 = str2double(answer{1});
    end
    close(h_fig)
    
    [m,n] = size(Cs_bis); %XXXXXXXXXXXXXXXXXXbetter to be based on bathy????????
    mean_Cs = mean(Cs_bis(~isnan(Cs_bis)));
    std_Cs = std(Cs_bis(~isnan(Cs_bis)));
    th1 = mean_Cs  + threshold_1 * std_Cs;
    
    for i=1:m
        for j=1:n
            if Cs_bis(i,j)>th1
               channels(i,j)=1;
            else 
               channels(i,j)=0;
            end
        end
    end

    channels(above_high_tide)=0; % No channels above high tide
    if casevars.isplot
        figure('Name','Network','Units','normalized','Tag','PlotFig');
        pcolor(channels)
        shading flat
        colorbar
        view([0 -90])
        set(gca,'clim',[-1 1])
        title('Extracted channel network')   
    end
end
%%
function results = get_number_channels(bathy_channels,casevars)
    %count the number of channels located at a specific distance from the coastal inlet
    % results(:,1) - Distance from the coastal inlet
    % results(:,2) - Number of channels
    % results(:,3) - Length of the semi-circle (m)
    
    results = [];
    cv = casevars;
    % CREATING CIRCLES AND COUNTING CHANNELS
    if cv.X_centre==0 || cv.Y_centre==0
       [cv,max_rad] = getXYcentre(bathy_channels,cv);
       if isempty(max_rad), return; end
    else
        [y_len, x_len] = size(bathy_channels);   
        my_rad = y_len-cv.Y_centre-1;
        mx_rad = min(cv.X_centre,x_len-cv.X_centre)-1;
        max_rad = min(mx_rad,my_rad); 
    end
    
    where_i=0;
    for radius = cv.min_rad:cv.rad_int:max_rad
        where_i=where_i+1;
        for x_i=0:2*radius+1
            x(x_i+1)=-radius+x_i;
        end
        for y_i=1:length(x)
            if y_i>radius+1
                y(y_i)=(radius^2-(x(y_i)-1)^2).^0.5;
            else
                y(y_i)=(radius^2-x(y_i)^2).^0.5;
            end
        end
        x_bathy=round(x+cv.X_centre);
        y_bathy=round(y+cv.Y_centre);
        n=1;
        profile=[];
        for i=1:length(x_bathy)
            profile(n)=bathy_channels(y_bathy(i),x_bathy(i));
            n=n+1;
            bathy_channels(y_bathy(i),x_bathy(i))=-1;
            if i<length(x_bathy)
               if y_bathy(i+1)-y_bathy(i)>1
                  dif=y_bathy(i+1)-y_bathy(i)-1;
                  half_dif=ceil(dif/2);
                  for dif_i=1:dif
                      if dif_i<=half_dif
                          profile(n)=bathy_channels(y_bathy(i)+dif_i,x_bathy(i));
                          n=n+1;
                          bathy_channels(y_bathy(i)+dif_i,x_bathy(i))=-1;
                      else
                          profile(n)=bathy_channels(y_bathy(i)+dif_i,x_bathy(i)+1);
                          n=n+1;
                          bathy_channels(y_bathy(i)+dif_i,x_bathy(i)+1)=-1;
                      end
                  end
               elseif y_bathy(i+1)-y_bathy(i)<-1
                  dif=y_bathy(i+1)-y_bathy(i)+1;
                  half_dif=ceil(dif/2);
                  for dif_i=-1:-1:dif
                      if dif_i>=half_dif
                         profile(n)=bathy_channels(y_bathy(i)+dif_i,x_bathy(i));
                         n=n+1;
                         bathy_channels(y_bathy(i)+dif_i,x_bathy(i))=-1;
                      else
                         profile(n)=bathy_channels(y_bathy(i)+dif_i,x_bathy(i)+1);
                         n=n+1;
                         bathy_channels(y_bathy(i)+dif_i,x_bathy(i)+1)=-1;
                      end
                  end

               end
            end       
        end
        count_channels=0;
        end_flag=1;
        for j=1:length(profile)
            if end_flag==1
               if profile(j)==1
                  count_channels=count_channels+1;
                  end_flag=0;
               end
            elseif end_flag==0
               if profile(j)==0
                  end_flag=1;
               end
            end
        end 

        results(where_i,1)=radius*cv.grid_size; % Distance from the coastal inlet
        results(where_i,2)=count_channels; % Number of channels
        results(where_i,3)=2*pi*radius*cv.grid_size/2; % Length of the semi-circle (m)
    end
    %
    if casevars.isplot
        figure('Name','Network','Units','normalized','Tag','PlotFig');
        pcolor(bathy_channels)
        shading flat
        colorbar
        view([0 -90])
        set(gca,'clim',[-1 1]);
        title('Radial circles on extracted grid')
    end
end
%%
function [cv,mxr] = getXYcentre(bathy_channels,cv)
    %Allow user to interactively adjust centre point
    [y_len, x_len] = size(bathy_channels);    
    if cv.X_centre==0
        cv.X_centre = round(x_len/2);
    end
    %
    if cv.Y_centre==0
        cv.Y_centre = round(y_len/2);
    end
    
    %plot base figure
    figtitle = sprintf('Centre selection');
    promptxt = 'Accept centre point definition';
    tag = 'PlotFig'; %used for collective deletes of a group
    butnames = {'Yes','No'};
    position = [0.2,0.4,0.4,0.4];
    [h_plt,h_but] = acceptfigure(figtitle,promptxt,tag,butnames,position);
    ax = axes(h_plt);
    pcolor(ax,bathy_channels)
    view([0 -90])
    shading flat
    hold on
    plot(ax,cv.X_centre,cv.Y_centre,'+r','MarkerSize',12)
    plot(ax,cv.X_centre,cv.Y_centre,'or','MarkerSize',10)
    xlabel('X-axis indices')
    ylabel('Y-axis indices')
    hold off
    ok = 0;
    while ok<1
        waitfor(h_but,'Tag');
        if ~ishandle(h_but) %this handles the user deleting figure window
            ok = 1;         %continue using default subdomain
        elseif strcmp(h_but.Tag,'No')
            %Get user to redfine subgrid
            promptxt = {'X centre','Y centre'};
            title = 'Define centre';
            defaults = string([cv.X_centre,cv.Y_centre]);
            asub = inputdlg(promptxt,title,1,defaults);
            newvals = cellfun(@str2double,asub)';
            cv.X_centre = newvals(1);
            cv.Y_centre = newvals(2);
            h_sd = findobj(ax,'Type','line');  %remove existing subdomain rectangle
            delete(h_sd);
            hold on
            plot(ax,cv.X_centre,cv.Y_centre,'+r','MarkerSize',12)
            plot(ax,cv.X_centre,cv.Y_centre,'or','MarkerSize',10)
            hold off
            h_but.Tag = '';
        else  %user accepted centre point - check not to near boundary
            my_rad = y_len-cv.Y_centre-1;
            mx_rad = min(cv.X_centre,x_len-cv.X_centre)-1;
            mxr = min(mx_rad,my_rad); 
            if mxr<cv.min_rad
                warndlg('Selected postion too close to edge to accept minimum radius')
                mxr = [];
                h_but.Tag = '';
            end
            delete(h_plt.Parent)
            ok = 1;
        end     
    end    
end
%%
function outable = plot_network_figures(bathy_0,channels,results,casevars)
    %generate the various plots to illustrate the network extracted and
    %resultant network properties
    figure('Name','Network properties','Units','normalized', ...                
                'Resize','on','HandleVisibility','on','Tag','PlotFig'); 

    s1 = subplot(2,2,1);
    pcolor(bathy_0)
    shading flat; 
    colorbar; 
    view([0 -90])
    set(s1,'clim',[-12 2]); 
    colormap jet 
    title('Input bathymetry')

    s2 = subplot(2,2,2);
    pcolor(channels)
    shading flat
    colorbar
    view([0 -90])
    set(s2,'clim',[-1 1])
    title('Extracted channel network')

    s3 = subplot(2,2,3);
    plot(s3,results(:,1),results(:,2),'o')
    xlabel ('Distance (m)')
    ylabel ('Number of channels')
    title(sprintf('Count of channels with\ndistance from entrance'))
    [maxcount,idx] = max(results(:,2));
    pktxt = sprintf('Peak count = %0.0f\nDistance = %0.0f',maxcount,results(idx,1));
    offset = s3.YLim(2)-(s3.YLim(2)-s3.YLim(1))*0.1;
    text(s3,s3.XLim(1)+1,offset,pktxt,'FontSize',8);

    s4 = subplot(2,2,4);
    AvDist = results(:,3)./results(:,2); %Av.drainage distance between channels (m)
    idinf = isinf(AvDist);
    AvDist(idinf) = [];  dist = results(~idinf,1);  %remove inf values from data
    nend = 3;                         %number of most distant points to omit
    plot(s4,dist(1:end-nend),AvDist(1:end-nend),'o')
    xlabel ('Distance (m)')
    ylabel (sprintf('Av. drainage distance\nbetween channels (m)'))  
    title(sprintf('Drainage distance with\ndistance from entrance'))
    mndist = mean(AvDist(1:end-nend));
    stdist = std(AvDist(1:end-nend));   
    [~,b,rsq,fitx,fity,fitxt] = regression_model(dist(1:end-nend),AvDist(1:end-nend),'linear');
    hold on
    xlimits = s4.XLim;
    plot(xlimits, mndist*[1 1],'--k');
    plot(xlimits, (mndist+stdist)*[1,1],':k');
    plot(xlimits, (mndist-stdist)*[1,1],':k');
    plot(s4,fitx,fity,'-.r');
    s4.XLim = xlimits;    %restore initial limits of data
    hold off
    offset = s4.YLim(1)+(s4.YLim(2)-s4.YLim(1))*0.1;
    text(s4,s4.XLim(1)+1,offset,fitxt,'FontSize',8);
    casedesc = sprintf('%s at %s',casevars.desc,casevars.time);  
    sgtitle(casedesc,'FontSize',12)
    
    %write summary statistics to the Commnad Window
    fprintf('Count = %0.0f, Peak distance = %0.0f, Drainage stats = %0.0f (%0.0f) [%0.3f(r^2=%0.2f)]\n',...
            maxcount,results(idx,1),mndist,stdist,b,rsq)
    outable = table(string(casedesc),maxcount,results(idx,1),mndist,stdist,b,rsq,...
               'VariableNames',{'CaseYear','PeakCount','PeakDistance',...
               'AvDrainDist','StdDrainDist','SlopeDrainDist','r2DrainDist'});
                    

    % figure;
    % plot(results(:,1),results(:,2),'--.')
    % xlabel ('Distance (m)')
    % ylabel ('Number of channels')
    % legend(casedesc{1})
    % title('Count of channels with distance from entrance')
end
%%
function Jd = geonet_diffusion(J,method,N,K,dt,sigma2) 
    % private function: diffusion (by Guy Gilboa):
    % Jd=diffusion(J,method,N,K)
    % Simulates N iterations of diffusion, parameters:
    %Input parameters
    % J =  source image (2D gray-level matrix) for diffusion
    % options - struct to define analysis parameters:
    %  method = 'lin':  Linear diffusion (constant c=1).
    %           'pm1': perona-malik, c=exp{-(|grad(J)|/K)^2} [PM90]
    %           'pm2': perona-malik, c=1/{1+(|grad(J)|/K)^2} [PM90]
    %           'tukey: c=0.5*((1-(grad(J)/K)^2)^2)*grad(J)
    %           'rmp': complex valued - ramp preserving [GSZ01]
    %  K    edge threshold parameter (lambda)
    %  N    number of iterations
    %  dt   time increment (0 < dt <= 0.25, default 0.2)
    %  sigma2 - if present, calculates gradients of diffusion coefficient
    %          convolved with gaussian of var sigma2 (Catte et al [CLMC92])J
    %Output parameters
    % Jd - image with selected diffusion smoothing applied

    %  This file is part of GeoNet.
    % 
    %     GeoNet is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    % 
    %     GeoNet is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    % 
    %     You should have received a copy of the GNU General Public License
    %     along with GeoNet.  If not, see <http://www.gnu.org/licenses/>.

    if ~exist('N','var')
       N=1;
    end
    if ~exist('K','var')
       K=1;
    end
    if ~exist('dt','var')
       dt=0.2;
    end
    if ~exist('sigma2','var')
       sigma2=0;
    end

    [Ny,Nx]=size(J); 

    if (nargin<3) 
       error('not enough arguments (at least 3 should be given)');
    end

    for i=1:N   
       % gaussian filter with kernel 5x5 (Catte et al)
       if (sigma2>0) 
          Jo = J;   % save J original
%           J=gauss(J,5,sigma2);  
          J = imgaussfilt(J,sigma2); %requires Matlab Image Processing Toolbox
       end

        % calculate gradient in all directions (N,S,E,W)
        In=[J(1,:); J(1:Ny-1,:)]-J;
        Is=[J(2:Ny,:); J(Ny,:)]-J;
        Ie=[J(:,2:Nx) J(:,Nx)]-J;
        Iw=[J(:,1) J(:,1:Nx-1)]-J;

        % calculate diffusion coefficients in all directions according to method
        switch method
            case 'lin'
                Cn=K; Cs=K; Ce=K; Cw=K;
            case 'pm1'
                Cn=exp(-(abs(In)/K).^2);
                Cs=exp(-(abs(Is)/K).^2);
                Ce=exp(-(abs(Ie)/K).^2);
                Cw=exp(-(abs(Iw)/K).^2);
            case 'pm2'
                Cn=1./(1+(abs(In)/K).^2);
                Cs=1./(1+(abs(Is)/K).^2);
                Ce=1./(1+(abs(Ie)/K).^2);
                Cw=1./(1+(abs(Iw)/K).^2);  
            case 'tukey'
                [lx, ly] = find(In < K);
                chr = sparse(lx, ly, 1, length(In), length(In(1,:)));
                Cn = .5 * ((1 - (In/K).^2).^2) .*chr;

                [lx, ly] = find(Is < K);
                chr = sparse(lx, ly, 1, length(Is), length(Is(1,:)));
                Cs = .5 * ((1 - (Is/K).^2).^2) .*chr;

                [lx, ly] = find(Ie < K);
                chr = sparse(lx, ly, 1, length(Ie), length(Ie(1,:)));
                Ce = .5 * ((1 - (Ie/K).^2).^2) .*chr;

                [lx, ly] = find(Iw < K);
                chr = sparse(lx, ly, 1, length(Iw), length(Iw(1,:)));
                Cw = .5 * ((1 - (Iw/K).^2).^2) .*chr;
            case 'rmp'
%                 k=K(1); theta=K(2); j=sqrt(-1);
                k=K; theta=K; j=sqrt(-1);
                Cn=exp(j*theta)./(1+(imag(In)/(k*theta)).^2);
                Cs=exp(j*theta)./(1+(imag(Is)/(k*theta)).^2);
                Ce=exp(j*theta)./(1+(imag(Ie)/(k*theta)).^2);
                Cw=exp(j*theta)./(1+(imag(Iw)/(k*theta)).^2);
            otherwise
                error(['Unknown method "' method '"']);
        end

       if (sigma2>0)   % calculate real gradiants (not smoothed)
        In=[Jo(1,:); Jo(1:Ny-1,:)]-Jo;
            Is=[Jo(2:Ny,:); Jo(Ny,:)]-Jo;
            Ie=[Jo(:,2:Nx) Jo(:,Nx)]-Jo;
          Iw=[Jo(:,1) Jo(:,1:Nx-1)]-Jo;
          J=Jo;
        end

       % Next Image J
       J=J+dt*(Cn.*In + Cs.*Is + Ce.*Ie + Cw.*Iw);
    end % for i

    Jd = J;
end
%%
function [xx,yy,qx,qy,mx,my] = geonet_qqplot(x,y,pvec)
    %QQPLOT Display an empirical quantile-quantile plot.
    %   QQPLOT(X) makes an empirical QQ-plot of the quantiles of
    %   the data set X versus the quantiles of a standard Normal distribution.
    %
    %   QQPLOT(X,Y) makes an empirical QQ-plot of the quantiles of
    %   the data set X versus the quantiles of the data set Y.
    %
    %   H = QQPLOT(X,Y,PVEC) allows you to specify the plotted quantiles in 
    %   the vector PVEC. H is a handle to the plotted lines. 
    %
    %   When both X and Y are input, the default quantiles are those of the 
    %   smaller data set.
    %
    %   The purpose of the quantile-quantile plot is to determine whether
    %   the sample in X is drawn from a Normal (i.e., Gaussian) distribution,
    %   or whether the samples in X and Y come from the same distribution
    %   type.  If the samples do come from the same distribution (same shape),
    %   even if one distribution is shifted and re-scaled from the other
    %   (different location and scale parameters), the plot will be linear.

    %   Copyright 1993-2000 The MathWorks, Inc. 
    %   $Revision: 2.12 $  $Date: 2000/05/26 18:53:31 $

    %  This file is part of GeoNet.
    % 
    %     GeoNet is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    % 
    %     GeoNet is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    % 
    %     You should have received a copy of the GNU General Public License
    %     along with GeoNet.  If not, see <http://www.gnu.org/licenses/>.

    if nargin == 1
       y  =  sort(x);
       [x,~]  = plotpos(y);
       x  = norminv(x); %requires Statistics and Machine Learning Toolbox.
       xx = x;
       yy = y;
    else
       n = -1;
       if nargin < 3
          nx = sum(~isnan(x));
          if (length(nx) > 1)
             nx = max(nx);
          end
          ny = sum(~isnan(y));
          if (length(ny) > 1)
             ny = max(ny);
          end
          n    = min(nx, ny);
          pvec = 100*((1:n) - 0.5) ./ n;
       end

       if (((size(x,1)==n) || (size(x,1)==1 && size(x,2)==n)) && ~any(isnan(x)))
          xx = sort(x);
       else
          xx = prctile(x,pvec);
       end
       if (((size(y,1)==n) || (size(y,1)==1 && size(y,2)==n)) && ~any(isnan(y)))
          yy = sort(y);
       else
          yy=prctile(y,pvec);
       end
    end

    q1x = prctile(x,25);
    q3x = prctile(x,75);
    q1y = prctile(y,25);
    q3y = prctile(y,75);
    qx = [q1x; q3x];
    qy = [q1y; q3y];


    dx = q3x - q1x;
    dy = q3y - q1y;
    slope = dy./dx;
    centerx = (q1x + q3x)/2;
    centery = (q1y + q3y)/2;
    maxx = max(x);
    minx = min(x);
    maxy = centery + slope.*(maxx - centerx);
    miny = centery - slope.*(centerx - minx);

    mx = [minx; maxx];
    my = [miny; maxy];

%     figure
%     hh = plot(xx,yy,'+',qx,qy,'-',mx,my,'-.');
%     if nargout == 1
%       h = hh;
%     end
% 
%     if nargin == 1
%        xlabel('Standard Normal Quantiles')
%        ylabel('Quantiles of Input Sample')
%        title ('QQ Plot of Sample Data versus Standard Normal')
%     else
%        xlabel('X Quantiles');
%        ylabel('Y Quantiles');
%     end

end
%%
function [pp,n] = plotpos(sx)
    %PLOTPOS Compute plotting positions for a probability plot
    %   PP = PLOTPOS(SX) compute the plotting positions for a probabilty
    %   plot of the columns of SX (or for SX itself if it is a vector).
    %   SX must be sorted before being passed into PLOTPOS.  The ith
    %   value of SX has plotting position (i-0.5)/n, where n is
    %   the number of rows of SX.  NaN values are removed before
    %   computing the plotting positions.
    %
    %   [PP,N] = PLOTPOS(SX) also returns N, the largest sample size
    %   among the columns of SX.  This N can be used to set axis limits.

    [n, m] = size(sx);
    if n == 1
       sx = sx';
       n = m;
       m = 1;
    end

    nvec = sum(~isnan(sx));
    pp = repmat((1:n)', 1, m);
    pp = (pp-.5) ./ repmat(nvec, n, 1);
    pp(isnan(sx)) = NaN;

    if (nargout > 1)
       n = max(nvec);  % sample size for setting axis limits
    end
end
%%
function C = curv(K,del)
    % Copyright (C) 2009  Passalacqua et al.

    %  This file is part of GeoNet.
    % 
    %     GeoNet is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    % 
    %     GeoNet is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    % 
    %     You should have received a copy of the GNU General Public License
    %     along with GeoNet.  If not, see <http://www.gnu.org/licenses/>.
    %K - ??
    %del - grid size

    [fx,fy] = gradient(K);
    fx=fx/del; 
    fy=fy/del; 
    fn = sqrt(fx.^2 + fy.^2);
    fx = fx ./ (fn + eps);   %0.000001
    fy = fy ./ (fn + eps);

    [fx1,fx2] = gradient(fx);
    [fy1,fy2] = gradient(fy);
    fx1=fx1/del; 
    fx2=fx2/del;
    fy1=fy1/del; 
    fy2=fy2/del; 
    C = fx1 + fy2;
end