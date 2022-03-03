%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Cases tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Cases*: deletes all Cases listed on the Cases tab but does not affect the model setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit the Case description.
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of Cases and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects Case(s) to be deleted from a list box of Cases and results are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.

%%
% <html>
% <table border=1><tr><td><u>Note</u>: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the <b>Project>Cases>Edit Data Set</b> option 
% to make a selection and then use the <i>Copy to Clipboard</i> button to paste 
% the selection to the clipboard.
% </td></tr></table>
% </html>

%% Setup > Input Data
% Load data from a file. To create a new instance (e.g. for a
% different location or data source) use Load. To add data to an
% existing data set, use Add.
%%
% * *Input Data > Gridded Data*: _Load/Add/Delete/Quality Control_ gridded data from file.
% * *Input Data > Timeseries*: _Load/Add/Delete/Quality Control_ timeseries data from file.
%%
% * _Load_: prompts user for file to be loaded. Once files have been read, user is prompted for a description (working title) for the data set. 
% * _Add_: prompts user for file to be loaded (only one file at a time can be added). Only files with the same dimensions as the existing data set can be used to Add data to a data record.
% * _Delete_: delete a grid from an existing Case table (ie a row).
% * _Quality Control_: data specific quality control to clean up data sets.
% Only available if defined for the specific data type.

%% Setup > Grid Tools
% * *Translate Grid*: interactively translate grid x-y
% coordinates.
% * *Rotate Grid*: interactively flip or rotate grid.   
% * *Re-Grid*: regrid a gridded dataset to match another grid or to user
% specified dimensions.
% * *Sub-Grid*: interactively define a subgrid and save grid as a new Case.               
% * *Combine Grids*: superimpose one grid on another based on maximum
% or minimum set of values.
% * *Add Surface*: add horizontal surface to an existing grid.
% * *To curvilinear*: map grid from cartesian to curvilinear coordinates. 
% * *From curvilinear*: map grid from curvilinear to cartesian coordinates.
% * *Difference Plot*: generate a plot of the difference between two grids.
% * *Export xyz Grid*: select a Case and export grid as xyz tuples.

%% Setup (other)
% * *Grid Parameters*: dialogue to set dimensions of default grid
% * *Run Parameters*: dialogue to define parameters for model skill estimates.
% * *Model Constants*: a number of constants are used in the model. Generally, the default values are appropriate but these can be adjusted and saved with the project if required.

%% Run
% * *Taylor Diagram*: generate a Taylor diagram using grids or timeseries
% relative to a reference grid or timeseries.
% * *Inlet Tools*:  tools to investigate grids which define inlets, 
% estuaries or tidal basins, including:
%%
% * >_Add Properties_: add form properties to a gridded data.
% * >_Delete Properties_: delete ALL property tables associated with a 
% selected gridded data set.
% * >_Tabulate Properties_: tabulate the gross properties for all timesteps in Selected Case.
% * >_Plot Properties_: plot selected gross properties for selected Cases as a function of time.
% * >_Plot Hypsometry_: plot surface area and volume hypsometry for one or more selected Cases and time intervals.
% * >_Plot X-Y property sets_: plot selected gross properties as X and Y to show co-variant change over time. Selected Cases must be the same length and are assumed to have contemporaneous points. 

%%
% * *User Tools*: options to generate additional properties, statistics and
% plots for gridded datasets, including analysis of networks (such as 
% dendritic channels) and rhythmic data (such as cusps).
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model.