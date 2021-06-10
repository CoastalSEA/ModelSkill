classdef MS_RunParams < muiPropertyUI     
%
%-------class help---------------------------------------------------------
% NAME
%   MS_RunParams.m
% PURPOSE
%   Class to handle ModelSkill run parameters
% USAGE
%   obj = MS_RunParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'X-axis definition ([x0 dx xN] or nint)',...
                          'Y-axis definition ([x0 dx xN] or nint)',...
                          'Maximimum correlation, R0',...
                          'Skill score exponent, n',...
                          'Skill score sub-sample window size, W',...
                          'Skill score averaging window (grids only)',...
                          'Skill score iteration option (0 or 1)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        Xaxis                %definition of the x-axis ([x0 dx xN] or nint)
        Yaxis                %definition of the y-axis ([x0 dx xN] or nint)
        maxcorr = 1          %maximimum correlation achievable (-)
        skillexponent = 1    %exponent to be used in skill score (-)
        skillwindow = 0      %number of points or grids to sub-sample over: +/-W (-)
        skillsubdomain       %used to average local skill over a sub domain
                             %subdomain defined as [x0,xN,y0,yN];
        skilliteration = 1   %flag to define iteration as
                             %true - iterates over every grid cell i=1:m-2W
                             %false - avoids overlaps and iterates over i=1:2W:m-2W
    end    

%%   
    methods (Access=protected)
        function obj = MS_RunParams(mobj)             % << Edit to classname
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'MS_RunParams';               % <<Edit to classname
            if isfield(mobj.Inputs,classname) && ...
                            isa(mobj.Inputs.(classname),classname)
                obj = mobj.Inputs.(classname);  
            else
                obj = MS_RunParams(mobj);             % << Edit to classname
            end
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            mobj.Inputs.(classname) = obj;
        end     
    end
%%        
        %add other functions to operate on properties as required   
end