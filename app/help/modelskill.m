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
% The ModelSkill App enables data to be loaded and then compared on a 
% Taylor diagram. This form of plot was originally proposed for the 
% comparison of model timeseries output (Taylor, 2001) and has subsequently
% been adapted for the comparison of morphological model outputs 
% (Bosboom and Reniers, 2014; Bosboom et al., 2014). In this implementation
% the Taylor approach is used but modified in line with the mehtod proposed by 
% Bosboom for the analysis of grid and mesh data. The options include the 
% ability to create and add points to a Taylor diagram, output the results 
% to the Clipboard, and estimate global and local skill scores.

%% Quick Startup
%  Quick Guide to setting up a new project
%  File>New 
%     To create a new project space.
% 	
%  Setup>Input Data>Gridded Data>Load
%     Load grids to be examined. Cases loaded are displayed on the Scenarios tab and the input files
%     used can be seen on the Inputs tab.
%  OR
%  Setup>Input Data>Timeseries>Load
%     Load timesereis data sets to be examined. Cases loaded are displayed on the Scenarios tab
%     and the input files used can be seen on the Inputs tab.	
% 
%  Setup>Grid Parameters
%     Define grid dimensions for use in when changing the gird dimensions using ReGrid
% 
%  Setup>run Parameters
%     Define the parameters to be used in the skill score - see Manual for details.
%     
%  File>Save as
%     Save model setup to a *.mat file.
% 	
%  Run>Taylor Diagram
%     Generates a Taylor diagram for selected Reference and Test data sets. Multiple data sets can
%     be added. Tabulates skill score and results can be copied to clipboard.
%
%  Run>Derive output
%     Manipulate data using external functions of Matlab functions to
%     created new datasets.
% 	
%  Plot>Plot Menu
%     Access standard plotting UI for 2d and 3d plotting options.
% 	
% Model runs are listed based on the user descriptions on the Cases tab.

%% Manual
% The <matlab:ms_open_manual manual> provides further details of setup and 
% configuration of the model. The files for the example use case can be found in
% the example folder <matlab:ms_example_folder here>. 

%% References
% Bosboom J and Reniers A, 2018, The Deceptive Simplicity of the Brier Skill Score, In: Handbook of Coastal and Ocean Engineering, Series, pp. 1639-1663.
%
% Bosboom J and Reniers A J H M, 2014, Scale-selective validation of morphodynamic models, 34th International Conference on Coastal Engineering, pp. 1911–1920, Seoul, South-Korea.
%
% Bosboom J, Reniers A J H M and Luijendijk A P, 2014, On the perception of morphodynamic model skill. Coastal Engineering, 94, 112-125, https://doi.org/10.1016/j.coastaleng.2014.08.008.
%
% Taylor K E, 2001, Summarizing multiple aspects of model performance in a single diagram. Journal of Geophysical Research - Atmospheres, 106 (D7), 7183-7192, 10.1029/2000JD900719.

%% See Also
% <matlab:doc('muitoolbox') muitoolbox>, <matlab:doc('dstoolbox') dstoolbox>.
	