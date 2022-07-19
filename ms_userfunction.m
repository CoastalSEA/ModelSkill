function hyps = ms_userfunction(mobj)
%
%-------function help------------------------------------------------------
% NAME
%   ms_userfunction.m
% PURPOSE
%   user function called from Inlet Tools menu
% USAGE
%   ms_userfunction(mobj)
% INPUTS
%   mobj - ModelUI instance
% OUTPUT
%   User defined
%   example returns hyps - a struct generated by gd_basin_volumes.m
% NOTES
%   The example code below calls gd_basin_hypsometry.m and gd_basin_volumes
%   to examine different ways of computing the volume between high and low
%   water (prism) in the basin.
% SEE ALSO
%   getInletTools. Called in ModelSkill
%
% Author: Ian Townend
% CoastalSEA (c) July 2022
%--------------------------------------------------------------------------
%

    %select case and timestep to analyse
    gridclasses = {'GD_ImportData'}; %add other classes if needed
    promptxt = {'Select Case to use:','Select timestep:'};
    [obj,~,irec] = selectCaseDatasetRow(mobj.Cases,[],...
                                                 gridclasses,promptxt,1);
    %get grid and water levels to use comput hypsometry
    grid = getGrid(obj,irec);
    [wl,histint] = gridWaterLevels(obj,true);                                         
    [~,hdst] = gd_basin_hypsometry(grid,wl,histint,0,true);
    %get summary volumes and plot
    hyps = gd_basin_volumes(hdst,true); %true to generate plots
    
    %check plots to compare zsurf/zvol from gd_basin_hypsometry (which is
    %the same as zsurf/zvol from gd_channel_hypsometry) with xzsurf and
    %xzvol derived values - use offset to make difference visible
    offset = 0;
    zsurf = sum(hyps.xzsurf,1);
    figure; plot(zsurf+offset,hyps.zcentre,hyps.zsurf,hyps.zcentre);
    hold on
    zvol = sum(hyps.xzvol,1);
    plot(zvol+offset,hyps.zcentre,hyps.zvol,hyps.zcentre);
    hold off
    ylabel('Elevation')
    xlabel('Volume or Surface Area')
    legend({'xzsurf','zsurf','xzvol','zvol'})
    subtitle("Compare zsurf/zvol from gd\_basin\_hypsometry with xzsurf and xzvol derived values")    
    
    %check prism derived from volumes with estimate using CSA as computed
    %in gd_section_properties
    for i=1:size(hyps.xzvol,1)
        zvolx = sum(hyps.xzvol(i:end,:),1);
        Vhw = interp1(hyps.zcentre,zvolx,wl.zhw(i),'linear','extrap');
        Vlw = interp1(hyps.zcentre,zvolx,wl.zlw(i),'linear','extrap');        
        Pr(1,i) = Vhw-Vlw;
    end  

    ps = gd_section_properties(grid,wl,hdst);
    delx = abs(grid.x(2)-grid.x(1));
    ixM = floor(grid.xM/(delx))+1;
    fprintf('prA = %.2e; prV = %.2e; vpr =  %.2e;\n',ps.PrA(ixM),ps.PrV(ixM),Pr(ixM))
    
    figure; 
    subplot(2,1,1)
    plot(hyps.xdist,ps.PrA); hold on, plot(hyps.xdist,Pr); hold off
    ylabel('Prism')
    legend({'Section properties prism','Basin volume prism'})
    prdiff = ps.PrA-Pr;
    subplot(2,1,2)
    plot(hyps.xdist,prdiff); 
    ylabel('Difference (Pr_Area-Pr_Volume)')
    xlabel('Distance')
    sgtitle('Comparison of along-channel tidal prism')
end