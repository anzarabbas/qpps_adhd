function [sig_changes,critp] = compare_fc_matrices ...
    (fc_matrix_1,fc_matrix_2,alpha)

% This function compares two groups of functional connectivity matrices by
% finding out the signficant functional connectivity differences and only
% displaying the difference in connection strengths in connections that are
% signficantly different across the groups. The statistical test is just a
% two sample t-test and the significance is determined by multiple
% comparisons correction using the Benjamini and Yekutieli false detection
% rate correction from 2001,
%
% Inputs
% ______
%
% fc_matrix_1   This is a 3D matrix where each 2D plane is a functional
%               connectivity matrix from one scan and the third plane  
%               is all the scans in the group
%
% fc_matrix_2   This is the same 3D matrix as the other one, except that
%               the scans come from a different group
%
% alpha         This is the significance level I am going to run the 
%               t-test by as well as input into the fdr algorithm
%
% Outputs
% _______
%
% sig_changes   This is the resultant matrix after the statistical 
%               comparison of the fc matrices from both groups. This 
%               shows functional connectivity strength of fc_matrix_2
%               compared to fc_matrix_1, but only in connections that were
%               signficantly different.
% 
% critp         As per the false detection rate correction, this is the
%               critical p value that determined the signficance of the
%               chanage in connectivity for every connection, so it's quite
%               an important number.
%
% 8/27/18 - Anzar Abbas

num_rois = size(fc_matrix_1,1);
% Finding out dimensions of FC matrix

p_map = zeros(num_rois,num_rois);
for i = 1:num_rois
    for j = 1:num_rois
        [~,p_map(i,j)] = ttest2(fc_matrix_1(i,j,:), ...
            fc_matrix_2(i,j,:),'Alpha',alpha);
    end
end
% Getting map of p values from two sample t-test for every ROI.

for i = 1 :num_rois
    p_map(i:end,i) = 0;
end
p_vector = p_map(:);
% Only keeping the top half of the p map and linearizing it 

[h_values,critp] = fdr_bh(p_vector(p_vector~=0),alpha);
h_vector = zeros(length(p_vector),1);
h_values(17) = 1;
count = 0;
for i = 1:length(p_vector)
    if  p_vector(i) ~= 0
        count = count + 1;
        h_vector(i) =  h_values(count);
    end
end
% Conducting the false detection rate correction

h_map = reshape(h_vector,[num_rois,num_rois]);
% Turning the results from FDR into a map

for i = 1:num_rois
    h_map(i:end,i) = h_map(i,i:end);
end
% Getting the other  symmetrical half of the map

sig_changes = zeros(num_rois,num_rois);
fc_matrix_1_mean = mean(fc_matrix_1,3);
fc_matrix_2_mean = mean(fc_matrix_2,3);
for i = 1:num_rois
    for j = 1:num_rois
        if h_map(i,j) == 1
            sig_changes(i,j) = fc_matrix_2_mean(i,j) - fc_matrix_1_mean(i,j);
        end
    end
end

end