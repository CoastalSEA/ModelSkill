%% ModelSkill Classes and Functions
% The ModelSkill App is built using the <matlab:doc('muitoolbox') muitoolbox>
% and a number of model specific classes.

%% ModelSkill classes
% * *ModelSkill*: defines the behaviour of the main UI.
% * *MS_RunParams*: UI to set pareameters used for the skill score
% calculation.

%% ModelSkill functions
% Functions to analyse univariate timeseries, or gridded data, using a Taylor diagram, 
% skill score and network analysis, along with toold to analyse the gross 
% properties (volume, area, width, etc) of inlets or channels 
%%
% * *example/demo_format*
% – used to set up import format and metadata description of timeseries variables. 
% * *getTaylorPlot* 
% - calls the taylor_plot function with either gridded or timeseries data.
% * *getInletTools* 
% - tools to extract and plot morphological properties of inlets and channels.
% * *getUserTools* 
% - user functions to do additional analysis on data loaded in ModelSkill. 
% * *ms_userfunction*
% - user function called from Inlet Tools menu. Example code compares
% different ways of computing basin tidal prism.
% * *network_count*
% - extract channels from a bathymetry and perform some network analysis.

%% Grid Classes
% Classes used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *GD_GridProps*: class inherits <matlab:doc('muipropertyui') muiPropertyUI> 
% abstract class, providing an interface to define the extent and intervals
% of a cartesian grid. 
% * *FGDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits *GDinterface*, which
% in turn inherits <matlab:doc('muidataset') muiDataSet> 
% to provide an extensive set of methods to handle datasets
% of various types (eg from models or imported files). 
% * *GD_ImportData*: class inherits <matlab:doc('fgdinterface') FGDinterface> abstract class (see above)
% to load xyz data from a file.

%% Grid Functions
% Functions used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *gd_ax_dir*
% - check direction of grid axes and reverse if descending, OR
% find grid orientation using ishead and direction of x-axis, OR
% check a grid axis direction by prompting user.
% * *gd_basin_hypsometry*
% - compute area and volume hypsometry from gridded elevation data.
% * *gd_basin_indices*
% - get the indices of the grid x-axis that fall within the basin or channel,
% when the mouth is offset from the grid origin. (NB: assumes basin/channel
% is aligned iwth the x-axis). Also returns the index of mouth position on 
% the x-axis.
% * *gd_basin_properties*
% - use the basin hypsometry from gd_basin_hypsometry to compute several 
% along-channel/x-axis morphological properties.
% * *gd_colormap*
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>) if not available.
% * *gd_digitisepoints*
% - creates figure to interactively digitise points on a grid and add
% elevations if required.
% * *gd_dimensions*
% - get the grid dimsnions for a grid struct (as used in GDinterface).
% * *gd_grid_line*
% - create a grid from scattered data points input as xyz tuples.
% * *gd_gross_properties*
% - compute the gross properties of a gridded bathymetry.
% * *gd_plan_form* 
% - compute planform variation along the x-axis at specified planar levels.
% * *gd_plotgrid*
% - create pcolor plot of gridded surface.
% * *gd_plotsections*
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
% * *gd_property_plots*
% - plots displayed on Proprety tab or stand-alone figure in Apps that use 
% GDinterface, such as ChannelForm and ModelSkill.
% * *gd_section_properties*
% - compute the width, cross-sectional area and prism along channel.
% * *gd_selectpoints*
% - accept figure to interactively select one or more points on a grid.
% * *gd_setpoint*
% - interactively select a point on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. for elevation).
% * *gd_startendpoints*
% - accept figure to interactively select start and end points on a grid.
% * *gd_subdomain*
% - accept figure to interactively select a subdomain of a grid.
% * *gd_property_plots* 
% - plots displayed on Proprety tab in ChannelForm model and on a figure 
% in ModelSkill.
% * *gd_xy2sn*
% - map grid from cartesian to curvilinear coordinates with option to return 
% the elevations on the source cartesian grid, or as a curvilinear grid.
% * *gd_sn2xy*
% - map grid from curvilinear to cartesian coordinates.
% * *getconvergencelength* 
% - least squares fit using fminsearch to % find the convergence length of 
% a channel from a distance-width xy data set.
% * *getsubgrid*
% - extract a subdomain from a grid and return the extracted grid and the 
% source grid indices of the bounding rectangle.
% * *a_star*
% - implements the A* search algorithm to find the shortest path given
% constraints (inaccessible cells) and a cost function (e.g. water depths).
% Author: Alex Ranaldi, 2022, https://github.com/alexranaldi/A_STAR
% * *InterX* 
% - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.
% * *xy2sn* 
% - Bart Vermeulen,2022, Cartesian to Curvilinear 
%   coordinate forward and backward transformation. 
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
% * *sn2xy* 
% - as above.

%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 

%% See Also 
% The ModelSkill <matlab:ms_open_manual manual>, which provides further details 
% of setup and configuration of the model.
