function [label_mask] = SubDivideImage(input_image, n_regions)
%% [label_mask] = SubDivideImage(input_image, n_regions)
%   Read in the image and ask the user to
if nargin == 0
    input_image = imread('test_pipe.tif');
    n_regions = input('How many regions? ');
elseif nargin < 2
    n_regions = input('How many regions? ');
end

imshow(input_image);
title('Input mask');

if n_regions == 0
    label_mask = zeros(size(input_image));
    return;
end
for r = 1:n_regions
    fprintf('Select regions for type %d\n',r);
    % make a 3D array of labels
    region(:,:,r) = r*roipoly;
end
nanregion = region;
nanregion(region==0) = nan;
label_mask = min(nanregion,[],3,'omitnan');
label_mask(isnan(label_mask)) = 0;
