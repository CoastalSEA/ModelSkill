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
% *NB*: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the *Project>Cases>Edit Data Set* option 
% to make a selection and then use the ‘Copy to Clipboard’ button to paste 
% the selection to the clipboard.

%% Setup
% * *Input Data*: import data from a file.
% * *Input Data > Gridded Data*: Load/Add/Delete/QC gridded data from file.
% * *Input Data > Timeseries*: Load/Add/Delete/QC timeseries data from file.
% * *Grid Parameters*: dialogue to define model run parameters.
% * *Grid Tools*: functions to manipulate a gridded dataset.
% * *Grid Tools > Translate Grid*: change the co-ordinates of the grid origin.
% * *Grid Tools > Rotate Grid*: flip or rotate the grid orientation.
% * *Grid Tools > Re-Grid*: Use defined grid dimensions or an existing
% grid to change the resolution of the grid (useful for doing differences).
% * *Grid Tools > Sub-Grid*: extract a sub-grid from an existing grid.
% * *Grid Tools > Combine Grids*: merge two grids based on the minimum or
% maximum values at each point in the grid.
% * *Grid Tools > Export xyz Grid*: export a grid as x,y,z tuples.
% * *Run Parameters*: dialogue to define parameters for model skill estimates.
% * *Model Constants*: a number of constants are used in the model. Generally, the default values are appropriate but these can be adjusted and saved with the project if required.

%% Run
% * *Taylor Diagram*: generate a Taylor diagram using grids or timeseries
% relative to a reference grid or timeseries.
% * *User Tools*: options to generate addtional properties or statistics
% for gridded datasets, including Gross properties (tables and plots),
% Hypsometry, Network analysis and a plot specific to rhythmic data (such as cusps).
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model.