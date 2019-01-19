
function [tc1, tc2] = get_roi_tc(roi1, func, roi2)

% This function calculates the timecourse of a region of interest (ROI) in
% the brain. Timecourse is defined as the mean signal at each timepoint.
%
% Inputs
% ______
%
% roi1      ROI mask. A 2D binary matrix of a region of interest in the
%           brain
%
% func      Functional timeseries. A 3D matrix in which the first two
%           dimensions are space and the third dimension is time. 
%
% roi2      Option input for a second ROI mask. A 2D binary matrix of a
%           region of interest in the brain unique from roi1.
%
% Outputs
% _______
%
% tc1       Timecourse of roi1. The mean signal of voxels in func from the
%           region defined in roi1 over time. This is a 1D vector.
%
% tc2       Timecourse of roi2 if inputted. The mean signal of voxels in
%           func from the region defined in roi2 over time. This is a 1D
%           vector.
%
% 7/28/16 - Anzar Abbas 

ts1 = zeros(size(func));
if nargin > 2
    ts2 = zeros(size(func));
end
% Predefining the matrices that will hold the timeseries for ROI1 and ROI2
% if specified.

t = size(func,3);
% Getting the length of the functional timeseries

for i = 1:t
    % For every timepoint
    
    ts1(:,:,i) = func(:,:,i) .* roi1;
    if nargin > 2
        ts2(:,:,i) = func(:,:,i) .* roi2;
    end
    % Calculating the timeseries from each ROI

end

tc1 = zeros(1,t);
if nargin > 2
    tc2 = zeros(1,t);
end
% Predefining the output timecourse vectors

for i = 1:t
    % For every timepoint
    
    instance = ts1(:,:,i);
    % Getting the image from one timepoint of the ROIs timeseries
    
    instance_sum = sum(instance(:));
    % Summing the signal from all the voxels in that timepoint 
    
    instance_length = length(find(instance(:)));
    % Finding out how many voxels are in there
    
    tc = instance_sum/instance_length;
    tc(isnan(tc)) = 0;
    tc1(i) = tc;
    % Calculating the mean of that timepoint manually
    % Also making sure there are no NaNs in there
    
    if nargin > 2
        instance = ts2(:,:,i);
        instance_sum = sum(instance(:));
        instance_length = length(find(instance(:)));
        tc = instance_sum/instance_length;
        tc(isnan(tc)) = 0;
        tc2(i) = tc;
    end
    % Doing the same for the second ROI if it was defined 
    
end

end
