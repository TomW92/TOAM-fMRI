function [ current ] = slice_time_correction( config )
%slice_time_correction sets up temporal scan
%
% Syntax: [ current ] = slice_time_correction( config )
%
% Inputs:
% config - config setup containing all information
%
% Outputs:
% current - spm slicetime output
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: config_setup.m
% Author: Mirko Bristle, Tom Willems
    
    current.temporal.st.scans = {config.current_session};
    current.temporal.st.nslices = config.params.nslices;
    current.temporal.st.tr = config.params.TR;
    current.temporal.st.ta = config.params.TA;
    current.temporal.st.so = config.params.slice_order;
    current.temporal.st.refslice = config.params.ref_slice;
    current.temporal.st.prefix = 'a';
end

