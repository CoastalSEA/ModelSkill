# ModelSkill
Matlab App for exploring the skill score between timeseries data sets or
gridded data sets.

## Licence
The code is provided as Open Source code (issued under a BSD 3-clause License).

## Requirements
ModelSkill is written in Matlab(TM) and requires v2016b, or later. In addition, ModelSkill requires the _dstoolbox_ and the _muitoolbox_.

## Background
The ModelSkill App enables data to be loaded and then compared on a Taylor diagram. This form of plot was originally proposed for the 
comparison of model timeseries output (Taylor, 2001) and has subsequently been adapted for the comparison of morphological model outputs (Bosboom and Reniers, 2014; Bosboom et al., 2014). In this implementation the Taylor approach is used but modified in line with the mehtod proposed by Bosboom for the analysis of grid and mesh data. The options include the ability to create and add points to a Taylor diagram, output the results to the Clipboard, and estimate global and local skill scores.

## References
Bosboom J and Reniers A, 2018, The Deceptive Simplicity of the Brier Skill Score, In: Handbook of Coastal and Ocean Engineering, Series, pp. 1639-1663.

Bosboom J and Reniers A J H M, 2014, Scale-selective validation of morphodynamic models, 34th International Conference on Coastal Engineering, pp. 1911–1920, Seoul, South-Korea.

Bosboom J, Reniers A J H M and Luijendijk A P, 2014, On the perception of morphodynamic model skill. Coastal Engineering, 94, 112-125, https://doi.org/10.1016/j.coastaleng.2014.08.008.

Taylor K E, 2001, Summarizing multiple aspects of model performance in a single diagram. Journal of Geophysical Research - Atmospheres, 106 (D7), 7183-7192, 10.1029/2000JD900719.

## Manual
The ModelSkill manual in the app/doc folder provides further details of setup and configuration of the model. The files for the example use case can be found in
the app/example folder. 

## See Also
The repositories for _dstoolbox_, _muitoolbox_ and _muiAppLIb_.