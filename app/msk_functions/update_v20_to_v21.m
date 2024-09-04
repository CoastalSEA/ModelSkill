function update_v20_to_v21(obj)
%
%-------header-------------------------------------------------------------
% NAME
%   update_v20_to_v21.m 
% PURPOSE
%   update saved models from v2.0 to v2.1.
% USAGE
%   update_v20_to_v21(obj)
% INPUTS
%   obj - instance of model
% RESULTS
%   saved model updated from v2.0 to v2.1.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
%   changes Form to Grid and adds Lt and Rv to Grid struct
%   assumes data to be updated is GD_ImportData class
%   To use from command line, open ModelSkill using:
%   >>ms = ModelSkill;     and then call
%   >>update_v20_to_v21(ms)
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2024
%--------------------------------------------------------------------------
%


%script to update from ModelSkill v2.0 to ModelSkill v2.1

msdata = obj.Cases.DataSets.GD_ImportData;

for i=1:length(msdata)
    grd = msdata(i);

    nrec = height(grd.Data.Form);
    
    grd.Data.Grid = grd.Data.Form;
    grd.Data = rmfield(grd.Data,'Form');
    
    grd.Data.Grid.UserData.Lt(nrec) = 0;
    grd.Data.Grid.UserData.Rv(nrec) = struct('Hr',[],'Wr',[],'Ar',[]);
    msdata(i) = grd;
    clear grd
end

obj.Cases.DataSets.GD_ImportData = msdata;