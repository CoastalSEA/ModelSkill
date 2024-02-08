function modelskill_update(obj,oldV,newV)
%
%-------header-------------------------------------------------------------
% NAME
%   modelskill_update.m 
% PURPOSE
%   update saved models to newer versions of ModelSkill
% USAGE
%   modelskill_update(oldV,newV) 
% INPUTS
%   obj - instance of model
%   oldV - old version number as a character string
%   newV - new version number as a character string
% RESULTS
%   saved model updated to new version. If this is called from ModelSkill this
%   will be the version that is being run at the time.
% NOTES
%   Called in muiModelUI.loadModel when old and new version numbers do not
%   match.
%
% Author: Ian Townend
% CoastalSEA (c) Feb 2024 
%--------------------------------------------------------------------------
%
    if strcmp(oldV,'2.0') && strcmp(newV,'2.1')
        update_v20_to_v21(obj);
    end
end