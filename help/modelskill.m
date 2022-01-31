%% ModelSkill
% Matlab App for exploring the skill score between timeseries data sets or
% gridded data sets.

%% Licence
% The code is provided as Open Source code (issued under a GNU General 
% Public License).

%% Requirements
% ModelUI is written in Matlab(TM) and requires v2016b, or later. In addition, 
% ModelUI requires both the <matlab:doc('dstoolbox') dstoolbox> and the 
% <matlab:doc('muitoolbox') muitoolbox>

%% Background

%% ModelSkill classes
% *ModelSkill* - defines the behaviour of the main UI.

%% ModelSkill functions
% Compare data sets (timeseries or gridded data) using Taylor diagram 
% and skill score.
%%	
%  Quick Guide to setting up a new project
%  File>New 
%     To create a new project space.
% 	
%  Setup>Input Data>Grid/Mesh Data
%     Load grids to be examined. Cases loaded are displayed on the Scenarios tab and the input files
%     used can be seen on the Inputs tab.
%  Setup>Input Data>Regrid Data
%     Allows grid to be interpolated onto the same grid used by other data sets or a user defined
%     grid.
%  OR
%  Setup>Input Data>Timeseries
%     Load timesereis data sets to be examined. Cases loaded are displayed on the Scenarios tab
%     and the input files used can be seen on the Inputs tab.	
% 
%  Setup>Run Parameters
%     Define grid dimensions for use in Regrid and parameters to be used in the skill score
%     see Manual for details.
% 	
%  File>Save as
%     Save model setup to a *.mat file.
% 	
%  Run>Taylor Diagram
%     Generates a Taylor diagram for selected Reference and Test data sets. Multiple data sets can
%     be added. Tabulates skill score and results can be copied to clipboard.
%  Run>User Tools
%     Access to a set of bespoke tools for analysing tidal inlets and estuaries to compute, tabulat
%     and plot gross properties, hypsometry and network properties.
%  Run>Derive output
%     Manipulate data using external functions of Matlab functions to
%     created new datasets.
% 	
%  Plot>Plot Menu
%     Access standard plotting UI for 2d and 3d plotting options.
% 	
% Scenarios are listed based on the user descriptions on the Cases tab.
%% Manual
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model. The files for the example use case can be found in
% the example folder <matlab:example_folder here>. 

%% See Also
% <matlab:doc('muitoolbox') muitoolbox>, <matlab:doc('dstoolbox') dstoolbox>.
	