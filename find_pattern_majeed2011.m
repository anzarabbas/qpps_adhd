
function [patt, corr] = find_pattern_majeed2011(func, mask, wl, stp)

% This function runs the pattern finding algorithm described in Majeed et
% al. 2011. It searches for a repeating event (pattern) within a defined
% region (mask) of a functional timeseries (func) of a certain length (wl)
% and starting from a certain timepoint (stp).
%
% Inputs
% ______
%
% func      Functional timeseries in 3D matrix. The first two dimensions
%           are space and the third dimension is time.
%
% mask      A 2D binary mask. This is the region in the brian where the
%           algorithm will be searching for a repeating event/pattern.
% 
% wl        Window length. This is a variable defining the length of the
%           event being searched for in units of TR.
%
% stp       Starting timepoint. This is the starting timepoint to be
%           entered into the algorithm.
%
% Outputs 
% _______
%
% patt      Pattern. This is the repeating event in the functional
%           timeseries according to the algorithm. It is a 3D matrix in
%           which the first two dimensions are space and the third
%           dimension is time (wl).
%
% corr      Sliding window correlation of patt with func. This is a 1D
%           vector depicting the Pearson correlation of patt with func over
%           the course of the entire scan.

corr_threshold = [0.1 0.2];
number_iterations_thresh1 = 3;
max_iterations = 20;
number_subjects = 1;
szI = size(func);
nt = szI(end)/number_subjects;
mask = (sum(abs(func), 3) > 0).*mask;
ind = mask(:) > 0;
I_msk = reshape(func, szI(1)*szI(2), szI(3));
I_msk = I_msk(ind, :);
szI_msk = size(I_msk);
patt = zeros(szI_msk(1), wl);
for k = 1:length(stp)
    patt = patt + I_msk(:, stp(k):...
        stp(k)+wl-1);
end
patt = patt(:) - mean(patt(:)); 
patt = patt / sqrt(patt'*patt);
corr = zeros(szI_msk(2)-wl+1, 1); 
QPPmaxima = zeros(size(corr));
for u = 1:number_subjects
    for k = 1:nt-wl+1
        im_chunk = I_msk(:,(u-1)*nt+k:(u-1)*nt+k+wl-1); 
        im_chunk = im_chunk(:) - mean(im_chunk(:));
        corr((u-1)*nt+k) = im_chunk'*patt/...
            sqrt(im_chunk'*im_chunk);
    end
    sub_t = corr( (u-1)*nt+1:(u-1)*nt+nt-wl+1 );
    sl = circshift(sub_t,1);
    sr = circshift(sub_t,-1);
    QPPmaxima((u-1)*nt+1:(u-1)*nt+nt-wl+1) = ((sub_t>sl)&...
        (sub_t>sr)).*sub_t;
    QPPmaxima((u-1)*nt+1) = 0; QPPmaxima((u-1)*nt+nt-wl+1) = 0; 
end
corr = corr(:);
t_old = corr;
t_old_old = corr;
t_old_old_old = corr;
p = 1;
while p<=max_iterations
    corr = smooth(corr);
    if(p <= number_iterations_thresh1)
        thresh = corr_threshold(1);
    else
        thresh = corr_threshold(2);
    end
    stp = find((QPPmaxima) > thresh);
    
    if (length(stp) == 1)
        break;
    end
    patt = zeros(szI_msk(1), wl);
    for k = 1:length(stp)
        patt = patt + I_msk(:, stp(k):...
            stp(k)+wl-1);
    end
    patt = patt(:) - mean(patt(:)); patt = ...
        patt / sqrt(patt'*patt);
    for u = 1:number_subjects
        for k = 1:nt-wl+1
            im_chunk = I_msk(:,(u-1)*nt+k:(u-1)*nt+k+wl-1); 
            im_chunk = im_chunk(:) - mean(im_chunk(:));
            corr( (u-1)*nt+k ) = im_chunk'*patt / ...
                sqrt(im_chunk'*im_chunk);
        end
        sub_t = corr( (u-1)*nt+1:(u-1)*nt+nt-wl+1 );
        sl = circshift(sub_t,1);
        sr = circshift(sub_t,-1);
        QPPmaxima((u-1)*nt+1:(u-1)*nt+nt-wl+1)=...
            ((sub_t>sl)&(sub_t>sr)).*sub_t;
        QPPmaxima((u-1)*nt+1) = 0; 
        QPPmaxima((u-1)*nt+nt-wl+1) = 0; 
    end
    corr = corr(:);
    if (ncc(t_old, corr) > .9999) || (ncc(t_old_old, corr)...
            > .9999)|| (ncc(t_old_old_old, corr) > .9999)
        break;
    end
    t_old_old_old = t_old_old;
    t_old_old = t_old;
    t_old = corr;
    p = p+1;
end
patt = zeros(szI(1), szI(2), wl);
for k = 1:length(stp)
    patt = patt + func(:, :, ...
        stp(k):stp(k)+wl-1);
end
patt = patt / length(stp);
function z = ncc(X, Y)
X = X(:) - mean(X(:));
Y = Y(:) - mean(Y(:));
if(norm(X) == 0 || norm(Y) == 0)
    z = nan;
    return;
end
z = (X' * Y) / norm(X)/norm(Y);
end
end
