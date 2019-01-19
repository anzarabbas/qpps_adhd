% ____________________________________ % 
%                                      % 
%            get_ac_networks           % 
% ____________________________________ %

function [c_mask, ac_mask, fc_map, min_max_value, max_min_value] = ...
    get_ac_networks(roi, func, thresh)

% This function takes an ROI, calculates its mean timecourse during a
% functional timeseries, and correlates the ROI timecourse with the
% timecourse of every voxel in the brain to output a static functional
% connectivity map of that ROI. It then uses that functional connectivity
% map to create a mask of the regions highly correlated with the ROI. It
% can also create a mask of the regions highly anticorrelated with the ROI.
% 
% Inputs
% ______
%
% roi       Region of interest mask. A 2D binary matrix of a region of
%           interest in the brain.
%
% func      Functional timeseries. A 3D matrix in which the first two
%           dimensions are space and the third dimension is time.
%
% thresh    The percentage of voxels in the brain that are to be added to
%           the c_mask and ac_mask matrices.
%
% Outputs 
% _______
%
% c_mask    Mask of areas highly correlated with ROI. This is a 2D binary
%           matrix acquired by separating the top (thresh)% of voxels in
%           fc_map
%
% ac_mask   Mask of areas anticorrelated with ROI. This is a 2D binary
%           matrix acquired by separating the bottom (thresh)% of voxels in
%           fc_map
%
% fc_map    Functional connectivity map of ROI in func. This is an 2D
%           matrix with the same spatial dimensions as func showing the
%           correlation between the mean timecourse of ROI and the
%           timecourse of every voxel in the brain.
%
% 12/18/16 - Anzar Abbas

roi_tc = get_roi_tc(roi,func);
% Getting the mean timecourse of the ROI

fc_map = zeros(size(roi));
% Predefining one of our output matries

for i = 1:size(func,1)
    for j = 1:size(func,2)
        if any(func(i,j,:))
        % Going through every brain voxel
        
        voxel_tc = func(i,j,:);
        % Getting the timecourse of the voxel
        
        x = corrcoef(roi_tc,voxel_tc);
        fc_map(i,j) = x(2);
        % Calculating the Pearson correlation between the voxel timecourse
        % and the ROI timecourse and adding it to fc_map
        
        end
    end
end

fc_map_brain = fc_map(fc_map~=0);
% Removing all non-brain voxels from the image and vectorizing it

min_max_value = prctile(fc_map_brain,(100-thresh));
% This is the minimum value in the maximum (thresh)% of voxels

max_min_value = prctile(fc_map_brain,thresh);
% This is the maximum value in the minimum (thresh)% of voxels

c_mask = 1 .* (fc_map > min_max_value);
% Assigning the top (thresh)% of voxels to c_mask

ac_mask = 1 .* (fc_map < max_min_value);
for i = 1:size(ac_mask,1)
    for j = 1:size(ac_mask,2)
        if fc_map(i,j) == 0
            ac_mask(i,j) = 0;
        end
    end
end
% Assigning the bottom (thresh)% of voxels to ac_mask

end
