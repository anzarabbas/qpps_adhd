
function func_regressed = regress_pattern(func, patt, corr)

% This function can regress a spatiotemporal pattern from the functional
% timeseries that it occurrs in. The regression is being done using a
% general linear model.
%
% Inputs
% ______
%
% func      Functional timeseries as a 3D matrix. The first two
%           dimensions are space and the third dimension is time.
%
% patt      The spatiotemporal pattern that I am trying to regress. This is
%           also a 3D matrix in which the first two dimensions are space
%           and the third dimension is time. What's different here compared
%           to func is that the time dimension is significantly shorter.
%
% corr      This is a sliding window correlation vector of patt with func.
%
% Outputs
% _______
%
% func_regressed    Functional timeseries as a 3D matrix after regression
%                   of pattern. The first two dimensions are space and the
%                   third dimension is time.
%
% Note that patt and corr are both outputs of the find_pattern_majeed2011
% function.
%
% 3/25/17 - Anzar Abbas

[x,y,t] = size(func);
% Getting the dimensions of the functional scan

func = reshape(func,[x*y,t]);
% Reshaping the functional scan so that space is now only one dimension.
% funct is now a 2D matrix in which the first dimension (rows) is voxels
% and the second dimension (columns) is time.

patt_length = size(patt,3);
% Finding out the number of timepoints in patt.

patt = reshape(patt,[x*y,patt_length]);
% Reshaping the patt so that space is now only in one dimension. patt is
% now a 2D matrix in which the first dimension (rows) is voxels and the
% second dimension (columns) is time.

corr = [ corr ; zeros(patt_length-1,1) ]';
% Making a new corr vector that will be the same length as the
% functional timeseries. Right now it is patt_length shorter than t.

func_regressed = zeros(size(func));
% Predefining the regressed version of the functional scan that will be the
% main output of this function 

for i = 1:(x*y)
    % For every brain voxel
    
    if any(func(i,:))
        % (Skipping all non-brain areas)
        
        regressor = conv(corr,patt(i,:),'valid');
        % Creating a regressor by convolving corr with the signal in this
        % particular voxel in patt.
        
        regressor = zscore(regressor);
        % Z-scoring the regressor
        
        % The following steps is the GLM regression
        z = (func(i,patt_length:end))';
        betas = (regressor*regressor')\regressor*z;
        z = z - regressor' * betas;
        
        func_regressed(i,patt_length:end) = z';
        % Adding the regressed voxel signal into the functional regressed
        % matrix 
        
    end
end

func_regressed = reshape(func_regressed,[x, y, t]);
% Putting the functional regressed scan back into the desired dimensions. 
        
end
