function data = readmatfile(filename,gridints)
%
%-------function help------------------------------------------------------
% NAME
%   readmatfile.m
% PURPOSE
%   read a mat z array from a file and reformat so that it can be loaded 
%   as a grid
% USAGE
%   data = readmatfile(filename,gridints)
% INPUTS
%   filename - name of file to be read
%   gridints = 1x2 cell array of characters defining x and y grid intervals
% OUTPUT
%   data - data set read from z array mat file dand formated as 3 cell
%          arrays with vectors of x, y and z to match format used in 
%          readinputfile
% NOTES
%   called in GD_ImportData as alternative option to text file import
%
% Author: Ian Townend
% CoastalSEA (c) Oct 2024
%--------------------------------------------------------------------------
%
    S = load(filename,'-mat');
    dataname = fieldnames(S);
    xint = str2double(gridints{1});
    yint = str2double(gridints{2});
    if isscalar(dataname)
        indata = S.(dataname{1});
        [m,n] = size(indata);
        [X,Y] = meshgrid(0:xint:(m-1)*xint,0:yint:(n-1)*yint);
        % figure
        % surf(X,Y,indata')
        data{1} = reshape(X,[],1);
        data{2} = reshape(Y,[],1);
        data{3} = reshape(indata,[],1);
    else
        warndlg('Multiple arrays in file. Data not loaded')
        data = [];
    end
end