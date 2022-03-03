%% ModelSkill Classes and Functions
% The ModelSkill App is built using the <matlab:doc('muitoolbox') muitoolbox>
% and a number of model specific classes.

%% ModelSkill classes
% * *ModelSkill*: defines the behaviour of the main UI.
% * *MS_RunParams*: UI to set pareameters used for the skill score
% calculation.

%% ModelSkill functions
% Compare data sets (timeseries or gridded data) using Taylor diagram 
% and skill score.
%%
% * *getTaylorPlot* - 
% - 
% * *getInletTools* - 

% * *getUserTools* - 
% - 
% * *network_count*
% - 



%% Grid Classes and Functions
% Classes and functions used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *GD_GridProps*: class inherits <matlab:doc('muipropertyui') muiPropertyUI> 
% abstract class, providing an interface to define the extent and intervals
% of a cartesian grid. 
% * *GDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits <matlab:doc('muidataset') muiDataSet> 
% and together they provide an extensive set of methods to handle datasets
% of various types (eg from models or imported files). 
% * *GD_ImportData*: class inherits <matlab:doc('gdinterface') GDinterface> abstract class (see above)
% to load xyz data from a file.
%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 

%% See Also 
% <matlab:doc('channelform_functions') ModelSkill functions> used in 
% ModelSkill App and the <matlab:open_manual manual>, which provides further details 
% of setup and configuration of the model.
