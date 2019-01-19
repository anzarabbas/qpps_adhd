
function combined_matrix = combine_fc_matrices(fc_matrix_1,fc_matrix_2)

% Since functional connectivity matrices are symmetric, half of them is
% a waste of space in a figure. Hence, I tend to put two functional
% connectivity matrices in one square to show off the most information with
% the least amount of space. This function does that for me, and all it
% needs are the two matrices that I want to combine.
%
% Inputs
% ______
%
% fc_matrix_1   is the first 2D functional connectivity matrix
% 
% fc_matrix_2   is the second 2D functional connectivity matrix
%
% Outputs
% _______
%
% combined_matrix   is the combination of the two functional connectivity
%                   matrices with fc_matrix_1 in the bottom left and
%                   fc_matrix_2 on the top right of the combined_matrix
%
% 8/27/18 - Anzar Abbas

num_rois = size(fc_matrix_1,1);
combined_matrix = zeros(num_rois,num_rois);
for i = 1:num_rois
    combined_matrix(i+1:end,i) = fc_matrix_1(i+1:end,i);
    combined_matrix(i,i+1:end) = fc_matrix_2(i,i+1:end);
end
