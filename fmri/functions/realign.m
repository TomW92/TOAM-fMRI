function [ current ] = realign( files )
%realign sets up the realignment with fieldmap or classic 
% Syntax: [ current ] = realign( config )
%
% Inputs:
% config - struct of session, subject file configuration
%
% Outputs:
% current - spm realign configuration
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: seupt_config.m
% Author: Mirko Bristle, Tom Willems

    current.spatial.realign.estwrite.data = files;
    current.spatial.realign.estwrite.eoptions.quality = 0.9;
    current.spatial.realign.estwrite.eoptions.sep = 4;
    current.spatial.realign.estwrite.eoptions.fwhm = 5;
    current.spatial.realign.estwrite.eoptions.rtm = 1;
    current.spatial.realign.estwrite.eoptions.interp = 2;
    current.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    current.spatial.realign.estwrite.eoptions.weight = '';
    
    %reslice
    current.spatial.realign.estwrite.roptions.which = [2 1];
    current.spatial.realign.estwrite.roptions.interp = 4;
    current.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    current.spatial.realign.estwrite.roptions.mask = 1;
    current.spatial.realign.estwrite.roptions.prefix = 'r';

