%% Requires CircStat toolbox from FileExchange 
% Homepage: https://philippberens.wordpress.com/code/circstats/
% Download: http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

clear;
close all;

DOCK_WINDOWS = false;

if DOCK_WINDOWS == true
    set(0,'DefaultFigureWindowStyle','docked');
else
    set(0,'DefaultFigureWindowStyle','normal');
end

plexus_image_file = uigetfile('.tif','Choose plexus image file'); 
plexus_thresh_image = imread(plexus_image_file);
% 
% %% Get data files, assume file basename is the same, +'_wss' for the WSS data.
[cell_data_file,root_dir] = uigetfile('*.csv','Choose polarity data file');
wss_data_file = strrep(cell_data_file,'.csv','_wss.csv'); 
xy_details = load(fullfile(root_dir,cell_data_file));

%% This allows a label mask to be used to identify datapoints according to a region
%  Initialised as zeros if no subdivision is required
labelMask = SubDivideImage(plexus_thresh_image);
%% *******
% Note: "csvread" doesn't like blank lines, "load" doesn't seem to mind.
flowData = load(fullfile(root_dir,wss_data_file));

% Filter for very large values (typical max around 40?)
% TODO: 
%   - This should be scaled to pixel physical size
%   - Maybe pass as an argument or option via GUI?
flow_data_to_exclude = find(hypot(flowData(:,1),flowData(:,2)) > 150); 

% If there are data points to exclude...
if sum(flow_data_to_exclude) > 0
    % ...Set to [0 0] which gets excluded later
    flowData(flow_data_to_exclude,:) = [0 0];
end
x0 = xy_details(:,1);
y0 = xy_details(:,2);
dx = xy_details(:,3)-xy_details(:,1);
dy = xy_details(:,4)-xy_details(:,2);
data = [x0,y0,dx,dy,flowData];

% Exclude cells with [0 0] flow as not being on the plexus
[cellOnPlexus,~] = find(~(flowData(:,1) == 0 & flowData(:,2) == 0));
data = data(cellOnPlexus,:);

for g = 1:size(data,1)
    % Round the index to find the label designation of the nearest pixel
    grp(g) = labelMask(round(data(g,2)),round(data(g,1)));
end

scalarProduct = dot([data(:,3),data(:,4)],[data(:,5),data(:,6)],2);

figure(1)
num_bins_scalar_product = sqrt(length(scalarProduct));
hist(scalarProduct,num_bins_scalar_product);
title('Histogram of scalar product');

figure(2)
lw = 1;
imshow(uint8(plexus_thresh_image)*55 + 200);
title('WSS and Cell Polarity data');
hold on;
quiver(data(:,1),data(:,2),data(:,3),data(:,4),'b');
quiver(data(:,1),data(:,2),data(:,5),data(:,6),'r');
legend({'Cell polarity','WSS'});
plot(data(:,1),data(:,2),'k.','MarkerSize',2);
axis image ij;


hold off;

%%
[cellAngle,cellLength] = cart2pol(data(:,3),data(:,4));
[flowAngle,flowLength] = cart2pol(data(:,5),data(:,6));

cellFlowAngle = cellAngle - flowAngle;
cellFlowAngleDeg = cellFlowAngle*180/pi;

figure(4)
number_of_bins = max(7,round(sqrt(length(cellFlowAngle))));
h = rose(cellFlowAngle,number_of_bins);
x = get(h,'Xdata');
y = get(h,'Ydata');

g=patch(x,y,'b');


%% Save MAT file to user-specified directory, containing entire workspace
timestamp = datetime('now','Format','yyyyMMdd-HHmmss');
[~,file_part,~] = fileparts(cell_data_file);
output_dir = uigetdir('root_dir','Choose directory for output files');
out_mat_file_name = strcat(file_part,'_AllData_',char(timestamp),'.mat');
save(fullfile(output_dir,out_mat_file_name));

