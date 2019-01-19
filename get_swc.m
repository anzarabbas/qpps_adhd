
function swc = get_swc(func,patt)

% This function calculates the sliding window correlation of a short
% timeseries (patt) with a long timeseries (patt)
%
% Inputs
% ______
% 
% func      Long functional timeseries. A 3D matrix in which the first two
%           dimensions are space and the third dimension is time. An
%           example of this would be a functional scan.
% 
% patt      Short functional timeseries. A 3D matrix in which the first two
%           dimensions are space and the third dimension is time. An
%           example of this would be a pattern that occurs within a longer
%           functional timeseries.
%
% Output
% ______
%
% swc       Sliding window correlation of patt with func
%
% 1/30/17 - Anzar Abbas

func_length = size(func,3);
% Getting the number of timepoints in func

patt_length = size(patt,3);
% Getting the number of timepoints in patt

swc = zeros(func_length-patt_length+1,1);
% Predefining the output vector

patt_row = patt(:);
% Turning this into a 1D vector 

for i = 1:length(swc)
    
    f = func(:,:,i:i+patt_length-1);
    f = f(:);
    % Getting the chunk of func that is of interest to us in this moment
    % and linearizing it into a 1D vector
    
    x = corrcoef(f,patt_row); x = x(2);
    % Getting the correlation between the functional chunk and the pattern
    
    swc(i) = x;
    % Filling in the correlation value in the output variable
    
end

end
