function [] = PlexusReconstruction_nonGuide()
% A package to reconstruct a vessel network for flow simulation
% Consists of 3 tabs:
%   First tab:
%       Reads in a mask from a segmented network
%       and extracts a skeleton and vessel diameter information
%       to generate a 3D reconstruction.
%       Takes extracted plexus data and runs HemeLB flow simulator
%   Second tab:
%       Allows the user to delineate cell polarity
%       (in this case 2 points per cell: nucleus + Golgi apparatus).
%   Third tab:
%       Allows statistical analysis of cell polarity vs computed
%       wall shear stress.
%       Allows subdivision of the plexus in regions of interest for
%       individual analysis.

% Brute force close all open figure windows, probably a more elegant way?
close all;

% Create main GUI window and register the tabs.
fullScreenSize = get(0,'Screensize');
UI.main_GUI = figure('position',[0 0 fullScreenSize(3) fullScreenSize(4)],...
              'menubar','figure',...
              'numbertitle','off',...
              'name','PolNet: polarity network analyser',...
              'resize','on',...
              'Position',[50 50 1400 800]);
UI.tabgp = uitabgroup(UI.main_GUI,'Position',[.005 .005 .99 .99]);
UI.tab1 = uitab(UI.tabgp,'Title','Generate Flow Model and Simulate');
UI.tab2 = uitab(UI.tabgp,'Title','Draw Polarity Vectors');
UI.tab3 = uitab(UI.tabgp,'Title','Analyse Polarity vs Wall Shear Stress');

% Generic button parameters. Width and height is relative to window size.
button_width = 0.13;
button_height = 0.05;
button_font_size = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First Tab: Generate Flow Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI.image_window_1 = axes('Parent',UI.tab1,'Position',[0.05 0.05 0.7 0.9]);
UI.image_window_2 = axes('Parent',UI.tab2,'Position',[0.05 0.05 0.7 0.9]);
UI.image_window_3 = axes('Parent',UI.tab3,'Position',[0.05 0.05 0.7 0.9]);

UI.wait_message = uicontrol('parent',UI.tab1, ...
                            'style','text','unit','normalized',...
                            'visible', 'off',...
                            'position',[0.8 0.8-12*button_height button_width 4*button_height],...
                            'ForegroundColor', 'red', 'FontSize',14,...
                            'string','There''s a time consuming operation running on the background. Please wait...');

UI.pushbutton_tab1(1) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8 button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Open File',...
                    'callback',{@open_callback_tab1,UI});
                
UI.pushbutton_tab1(2) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-1*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Skeletonise',...
                    'callback',{@skeletonise_callback,UI});

UI.pushbutton_tab1(3) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-2*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Reconstruct surface',...
                    'callback',{@reconstruct_callback,UI});

UI.pushbutton_tab1(4) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-3*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Set up flow simulation',...
                    'callback',{@setup_HemeLB_callback,UI});

UI.pushbutton_tab1(5) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-4*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Run flow simulation',...
                    'callback',{@run_HemeLB_callback,UI});

UI.radio_button_tab1 = uicontrol('Parent',UI.tab1,...
                                 'Style','checkbox','unit','normalized',...
                                 'String','Run previously set up model',...
                                 'fontsize',button_font_size,...
                                 'Value',0,...
                                 'Position',[0.8 0.8-5*button_height 2*button_width button_height],...
                                 'callback', {@use_previous_model_setup_checkbox_callback, UI});
use_previous_model_setup_checkbox_callback(UI.radio_button_tab1, [], UI) % Force callback once in case user doesn't interact with control

UI.pushbutton_tab1(6) = uicontrol('Parent',UI.tab1,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-7*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Quit',...
                    'callback',{@quit_callback,UI});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second Tab: Draw Polarity Vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI.pushbutton_tab2(1) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8 button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Open File',...
                    'callback',{@open_callback_tab2,UI});

UI.pushbutton_tab2(2) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-5*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Quit',...
                    'callback',{@quit_callback,UI});

UI.pushbutton_tab2(3) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-2*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Start Selecting Cells',...
                    'Tag','start_selecting_cells',...%                     'Value',1.0,...
                    'callback',{@start_selecting_cells_callback,UI});
                
UI.pushbutton_tab2(4) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-1*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Load position data',...
                    'callback',{@load_position_data_callback,UI});
                
UI.pushbutton_tab2(5) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-3*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Save position data',...
                    'callback',{@save_position_data_callback,UI});
                
UI.pushbutton_tab2(6) = uicontrol('Parent',UI.tab2,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-4*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Delete Last Arrow',...
                    'callback',{@delete_last_vector_callback,UI}); 
   
UI.crosshair_text_tab2 = uicontrol('Parent',UI.tab2,...
                'Style','text','unit','normalized',...
                'fontsize',button_font_size,...
                'String','Crosshair colour',...
                'Position',[0.8 0.8-7*button_height 0.4*button_width button_height]);

UI.dropdown_tab2   = uicontrol('Parent',UI.tab2,...
                    'style','popup','unit','normalized',...
                    'position',[0.8+0.4*button_width 0.8-7*button_height 0.6*button_width button_height],...
                    'fontsize',button_font_size-2,...
                    'string',{'Magenta','White','Red','Green','Blue','Cyan','Black','Yellow'},...
                    'callback',{@dropdown_callback,UI}); 
 % Force dropdown callback once in case user doesn't interact with control
 dropdown_callback(UI.dropdown_tab2, [], UI)
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Third Tab: Analyse Polarity vs Wall Shear Stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UI.pushbutton_tab3(1) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8 button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Get flow info at cell nuclei',...
                    'callback',{@get_wss_nuclei,UI});

UI.dropdown_tab3 = uicontrol('Parent',UI.tab3,...
                             'style','popup','unit','normalized',...
                             'position',[0.8 0.8-button_height button_width button_height],...
                             'fontsize',button_font_size,...
                             'string',{'Wall shear stress','Velocity','Pressure'},...
                             'callback',{@flow_variable_dropdown_callback,UI});
% Force callback once in case user doesn't interact with control
flow_variable_dropdown_callback(UI.dropdown_tab3, [], UI)

UI.radio_button_tab3 = uicontrol('Parent',UI.tab3,...
                                 'Style','checkbox','unit','normalized',...
                                 'String','Use last simulation results',...
                                 'fontsize',button_font_size,...
                                 'Value',1,...
                                 'Position',[0.8 0.8-2*button_height 2*button_width button_height],...
                                 'callback', {@use_last_results_checkbox_callback, UI});
% Force callback once in case user doesn't interact with control
use_last_results_checkbox_callback(UI.radio_button_tab3, [], UI)

UI.pushbutton_tab3(2) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-4*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Display mask and vectors',...
                    'callback',{@show_vectors_callback,UI});
                
UI.pushbutton_tab3(3) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-5*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Subdivide into regions',...
                    'callback',{@subdivide_regions_callback,UI});

UI.pushbutton_tab3(4) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-6*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Analyse Regions',...
                    'callback',{@analyse_regions_callback,UI});

UI.pushbutton_tab3(5) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-7*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Save Data',...
                    'callback',{@save_analysis_data_callback,UI});

UI.pushbutton_tab3(6) = uicontrol('Parent',UI.tab3,...
                    'style','pushbutton','unit','normalized',...
                    'position',[0.8 0.8-8*button_height button_width button_height],...
                    'fontsize',button_font_size,...
                    'string','Quit',...
                    'callback',{@quit_callback,UI});
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% END OF GUI BUILDING FUNCTION             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

% CONSTANTS (the MATLAB way)
function [MAX_CELLS] = MAX_NUM_CELLS
    MAX_CELLS = 5000;

function [FREQUENCY] = AUTOSAVE_FREQUENCY
    FREQUENCY = 5;

function [MESSAGE] = ERROR_MESSAGE
    MESSAGE = 'Please check the PolNet terminal for error messages and refer to the Troubleshooting section of the paper.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% Callbacks                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

%% Callbacks for first tab
%--------------------------------------------------------------------------
function [] = quit_callback(hObject,eventdata,handles)
% Function used by all tabs
user_response = questdlg('Are you sure you want to close the program?');
switch user_response
case 'No'
    % take no action
case 'Yes'
    % Prepare to close GUI application window
    % Autosave state?
    delete(handles.main_GUI)
    % Close all the extra figures (TODO: make this specific)
    close all
end

%--------------------------------------------------------------------------
function [] = open_callback_tab1(hObject, eventdata, handles)
% Check valid file format
validFileTypes = {'*.tif*';'*.jpg'};
[inputFile,inputDir] = uigetfile(validFileTypes);
[~,filePrefix,~] = fileparts(inputFile);
set(handles.main_GUI,'currentaxes',handles.image_window_1);
fileToOpen = fullfile(inputDir,inputFile);
input_img = imread(fileToOpen);
% Convert to logical
input_img = im2bw(input_img);

imshow(input_img);
title(filePrefix,'Interpreter','none');
setappdata(handles.main_GUI,'plexus_mask_file',fileToOpen);
setappdata(handles.main_GUI,'plexus_mask_image',input_img);

%--------------------------------------------------------------------------
function [] = skeletonise_callback(hObject, eventdata,handles)
mask_img = getappdata(handles.main_GUI,'plexus_mask_image');
mask_file = getappdata(handles.main_GUI,'plexus_mask_file');

dialogOutput = inputdlg('Please provide pixel/micron ratio:');
if (isempty(dialogOutput))
    % User cancelled previous dialog
    return
end
pixelsPerUm = str2double(dialogOutput{1});
if isnan(pixelsPerUm)
    uiwait(msgbox('Wrong input. Try again.'));
    return;
end
set(handles.wait_message,'visible','on');
[skeleton_img, skeleton_file] = SkeletonizeTiffPlexus(pixelsPerUm, mask_file, mask_img);

setappdata(handles.main_GUI,'skeleton_file',skeleton_file);
% Overlay skeleton onto mask image in window
% dilate to >1 pixel for viewing, otherwise can't see on large images
dilate_radius = round(max(size(mask_img))/(2*1024)); % Try to make it visible in a 1024x1024 image
% Use De Morgan's law to view skeleton and mask at same time
% Arrange view as white skeleton on black mask, regardless of input
imgMax = max(mask_img(:));
corner(1) = mask_img(1,1);
corner(2) = mask_img(1,end);
corner(3) = mask_img(end,1);
corner(4) = mask_img(end,end);
if mean(corner) > imgMax/2 % It has white background
    skeleton_view = mask_img | imdilate(skeleton_img,strel('disk',dilate_radius));
else
    skeleton_view = ~mask_img | imdilate(skeleton_img,strel('disk',dilate_radius));
end
set(handles.main_GUI,'currentaxes',handles.image_window_1);
imshow(skeleton_view);
set(handles.wait_message,'visible','off');


%% Callbacks for second tab
%--------------------------------------------------------------------------
function [] = open_callback_tab2(hObject, eventdata, handles)
% Check valid file format
validFileTypes = {'*.tif*';'*.lsm'};
[inputFile,inputDir] = uigetfile(validFileTypes);
[~,filePrefix,~] = fileparts(inputFile);
set(handles.main_GUI,'currentaxes',handles.image_window_2);
imshow(fullfile(inputDir,inputFile));

title(filePrefix,'Interpreter','none'); % Stop it parsing as LaTeX for title

setappdata(handles.tab2,'nuclei_golgi_file_dirname',inputDir);

%--------------------------------------------------------------------------
function start_selecting_cells_callback(hObject, eventdata, handles)
% If we haven't collected or loaded any points yet
if ~isappdata(handles.image_window_2,'points_1')
    % Initialise empty arrays
    points_1 = cell(MAX_NUM_CELLS, 1);
    points_2 = cell(MAX_NUM_CELLS, 1);
    num_points = 0;
else 
    points_1 = getappdata(handles.image_window_2,'points_1');
    points_2 = getappdata(handles.image_window_2,'points_2');
    num_points = getappdata(handles.image_window_2,'num_points');
end

setappdata(gca,'num_points',num_points);

% Keep looping until user clicks any button other than left click
while 1
    num_points = num_points + 1;
    setappdata(gca,'num_points',num_points);

    [x1,y1,button1] = ginput_pointer(1);

    % If anything other than a left click, exit the loop without taking
    % point
    if button1 ~= 1 
        points_1{num_points} = [];
        points_2{num_points} = [];
        num_points = num_points - 1; % decrement counter
        setappdata(hObject,'string','Start Selecting Cells');
        setappdata(gca,'num_points',num_points);
        break;
    end
    points_1{num_points} = [x1,y1];
    ShowPointsPlot(points_1,'nucleus');

    [x2,y2,button2] = ginput_pointer(1);
    % If anything other than a left click, exit the loop without recording
    % point
    if button2 ~= 1
        points_1{num_points} = [];
        points_2{num_points} = [];
        num_points = num_points - 1; % decrement counter
        setappdata(hObject,'string','Start Selecting Cells');
        setappdata(gca,'num_points',num_points);
        % Need to replot to get rid of dangling nucleus point
        ShowPointsPlot(points_1,'nucleus'); 
        break;
    end
    points_2{num_points} = [x2,y2];
    ShowPointsPlot(points_2,'golgi');

    setappdata(handles.image_window_2,'points_1',points_1);
    setappdata(handles.image_window_2,'points_2',points_2);
    setappdata(handles.image_window_2,'num_points',num_points);

    ShowQuiverPlot

    if mod(num_points, AUTOSAVE_FREQUENCY) == 0
        save_position_data(handles, 'autosave_polarity_vectors.csv', getappdata(handles.tab2,'nuclei_golgi_file_dirname'));
    end
end

%--------------------------------------------------------------------------
function dropdown_callback(hObject, eventdata, handles)
colour_idx = get(hObject,'Value');
if isempty(colour_idx)
    colour_idx = 1;
end
current_axes = gca;
switch colour_idx
    case 1
        setappdata(current_axes,'crosshair_colour','m');
    case 2
        setappdata(current_axes,'crosshair_colour','w');
    case 3
        setappdata(current_axes,'crosshair_colour','r');
    case 4
        setappdata(current_axes,'crosshair_colour','g');
    case 5
        setappdata(current_axes,'crosshair_colour','b');
    case 6
        setappdata(current_axes,'crosshair_colour','c');
    case 7
        setappdata(current_axes,'crosshair_colour','k');
    case 8
        setappdata(current_axes,'crosshair_colour','y');        
    otherwise
        disp('Unknown colour') % Set to magenta 
        setappdata(current_axes,'crosshair_colour','m');
end

%--------------------------------------------------------------------------
function ShowQuiverPlot(~)
currentAxes = gca;
points_1 = getappdata(currentAxes,'points_1');
points_2 = getappdata(currentAxes,'points_2');
nCells = getappdata(currentAxes,'num_points');

% Quiver works in X,Y, not row,col
nucleus = zeros(nCells,2);
golgi = zeros(nCells,2);

for c = 1:nCells
    nucleus(c,:) = points_1{c};
    golgi(c,:) = points_2{c};
end
 
dxdy = golgi-nucleus;

% We don't want the vectors to overlap with the impoints, otherwise the
% latter become tricky to select with the mouse. Draw a vector that is 10%
% shorter than the actual distance between nucleus and golgi.
vectorOrigin = nucleus + 0.05*dxdy;
vectorLength = 0.9*dxdy;

% delete old quiver plot (if it exists) and update, in case any points were moved
delete(getappdata(currentAxes,'quiver'));
hold on;
newQuiver = quiver(currentAxes, vectorOrigin(:,1),vectorOrigin(:,2),vectorLength(:,1),vectorLength(:,2),0,'color','r');
setappdata(currentAxes,'quiver',newQuiver);
    
%--------------------------------------------------------------------------
function ShowPointsPlot(data_points,data_type)
num_points = getappdata(gca,'num_points');
if num_points == 0
    nucleus_plot = getappdata(gca,'nucleus_plot');
    golgi_plot = getappdata(gca,'golgi_plot');
    delete(nucleus_plot);
    delete(golgi_plot);
    return;
end
if strcmp(data_type,'nucleus')
    nucleusData = cat(1,data_points{1:num_points});
    if isappdata(gca,'nucleus_plot')
        delete(getappdata(gca,'nucleus_plot'));
    end
    hold on;
    nucleus_plot = plot(nucleusData(:,1),nucleusData(:,2),'b.','MarkerSize',30);
    setappdata(gca,'nucleus_plot',nucleus_plot);
elseif strcmp(data_type,'golgi')
    golgiData = cat(1,data_points{1:num_points});
    if isappdata(gca,'golgi_plot')
        delete(getappdata(gca,'golgi_plot'));
    end
    hold on;
    golgi_plot   = plot(golgiData(:,1),golgiData(:,2),'g.','MarkerSize',30);
    setappdata(gca,'golgi_plot',golgi_plot);
else
    disp('Unknown data type, expecting nucleus or golgi');
end

%--------------------------------------------------------------------------
function load_position_data_callback(hObject, eventdata, handles)
nuclei_golgi_file_dirname = getappdata(handles.tab2,'nuclei_golgi_file_dirname');
[dataFile,dataDir] = uigetfile([nuclei_golgi_file_dirname '*.csv']);
positionData = load(fullfile(dataDir,dataFile));
[nRows,nCols] = size(positionData);
if nCols ~= 4
    disp([dataFile ' has wrong file format: 4 columns expected.'])
    % Invalid file
    return 
end

hNucleus = cell(max(nRows, MAX_NUM_CELLS), 1);
hGolgi = cell(max(nRows, MAX_NUM_CELLS), 1);
for row = 1:nRows
    hNucleus{row} = [positionData(row,1),positionData(row,2)];
    hGolgi{row} = [positionData(row,3),positionData(row,4)];
end

setappdata(handles.image_window_2,'points_1',hNucleus);
setappdata(handles.image_window_2,'points_2',hGolgi);
setappdata(handles.image_window_2,'num_points',nRows);

ShowPointsPlot(hNucleus,'nucleus');
ShowPointsPlot(hGolgi,'golgi');
ShowQuiverPlot

%--------------------------------------------------------------------------
function save_position_data_callback(hObject, eventdata, handles)
nuclei_golgi_file_dirname = getappdata(handles.tab2,'nuclei_golgi_file_dirname')
[dataFile,dataDir] = uiputfile([nuclei_golgi_file_dirname '*.csv']);
save_position_data(handles, dataFile, dataDir)

%--------------------------------------------------------------------------
function save_position_data(handles, dataFile, dataDir)
points_1 = getappdata(handles.image_window_2,'points_1');
points_2 = getappdata(handles.image_window_2,'points_2');
nPoints = getappdata(handles.image_window_2,'num_points');

if nPoints == 0
    disp('There are no points to save!')
    return
end

data = zeros(nPoints,4);

for p = 1:nPoints
    data(p,:) = [points_1{p},points_2{p}];
end
dlmwrite(fullfile(dataDir,dataFile),data,'precision','%.2f')

%--------------------------------------------------------------------------
function [] = delete_last_vector_callback(hObject, eventdata, handles)
current_axes = gca;
num_points = getappdata(current_axes,'num_points');
if num_points == 0
    % There are no points to delete!
    disp('There are no points to delete!');
elseif num_points == 1
    % We're trying to delete the last vector
    % TODO: Grey out the button after deleting final vector
    num_points = 0;
    setappdata(current_axes,'num_points',num_points);
    ShowPointsPlot([],'nucleus');
    ShowPointsPlot([],'golgi');
    ShowQuiverPlot;
elseif num_points > 1
    points_1 = getappdata(current_axes,'points_1');
    points_2 = getappdata(current_axes,'points_2');
    
    num_points = getappdata(current_axes,'num_points');
    points_1{num_points} = [];
    points_2{num_points} = [];
    
    setappdata(current_axes,'points_1',points_1);
    setappdata(current_axes,'points_2',points_2);
    
    num_points = num_points - 1;
    setappdata(current_axes, 'num_points',num_points);
    ShowPointsPlot(points_1,'nucleus');
    ShowPointsPlot(points_2,'golgi');
    ShowQuiverPlot
end

%--------------------------------------------------------------------------
function [] = reconstruct_callback(hObject, eventdata,handles)
set(handles.wait_message,'visible','on');
skeleton_file = getappdata(handles.main_GUI,'skeleton_file');
[voxel_size, time_step, num_time_steps, stl_file_name] = ReconstructSurfaceFromSkeleton(skeleton_file);
pr2_file_name = generate_pr2_file(voxel_size, time_step, num_time_steps, stl_file_name);
setappdata(handles.main_GUI,'pr2_file_name', pr2_file_name);
set(handles.wait_message,'visible','off');

%--------------------------------------------------------------------------
function [pr2_file_name] = generate_pr2_file(voxel_size, time_step, num_time_steps, stl_file_name)
[stl_file_path, stl_file_basename, ~] = fileparts(stl_file_name);
pr2_file_name = fullfile(stl_file_path, [stl_file_basename '.pr2']);
file_id = fopen(pr2_file_name, 'w');
minimum_pr2_file = {'DurationSeconds', num2str(time_step*num_time_steps);
                    'Iolets', '[]';
                    'OutputGeometryFile', fullfile(stl_file_path, [stl_file_basename '.gmy']);
                    'OutputXmlFile', fullfile(stl_file_path, [stl_file_basename '.xml']);
                    'SeedPoint', '{x: .nan, y: .nan, z: .nan}';
                    'StlFile', stl_file_name;
                    'StlFileUnitId', '2';
                    'TimeStepSeconds', num2str(time_step);
                    'VoxelSize', num2str(voxel_size)};
[num_file_lines, ~] = size(minimum_pr2_file);
for line = 1:num_file_lines
    fprintf(file_id, '%s: %s\n', minimum_pr2_file{line,:});
end
fclose(file_id);

%--------------------------------------------------------------------------
function [] = setup_HemeLB_callback(hObject,eventdata,handles)
if ispc
    errordlg('I don''t know how to run HemeLB''s setup tool under Windows.')
else
    set(handles.wait_message,'visible','on');
    pause(0.1); % short pause ensures wait message is displayed before calling system()
    pr2_file_name = getappdata(handles.main_GUI,'pr2_file_name');
    command = ['hemelb-setup --pro ' pr2_file_name];
    ret = system(command);
    if ret ~= 0
        title = 'Flow simulation setup failed';
        uiwait(msgbox(ERROR_MESSAGE, title));
        set(handles.wait_message,'visible','off');
        return;
    end

    [pr2_file_path, pr2_file_basename, ~] = fileparts(pr2_file_name);
    xml_file_name = fullfile(pr2_file_path, [pr2_file_basename '.xml']);
    command = ['python ModifyXMLFile.py ' xml_file_name];
    ret = system(command);
    if ret ~= 0
        title = 'Flow simulation setup failed';
        uiwait(msgbox('Cannot modify .xml input file. Simulation may produce incomplete results.', title));
        set(handles.wait_message,'visible','off');
        return;
    end
    set(handles.wait_message,'visible','off');
end

%--------------------------------------------------------------------------
function [] = run_HemeLB_callback(hObject,eventdata,handles)
if ispc
    errordlg('I don''t know how to run HemeLB under Windows.')
else
    if getappdata(handles.tab1, 'use_previous_model_setup')
        [xml_file_name, xml_file_directory] = uigetfile('*.xml');
        xml_file_full_path = fullfile(xml_file_directory, xml_file_name);
    else
        pr2_file_name = getappdata(handles.main_GUI,'pr2_file_name');
        [xml_file_directory, xml_file_basename, ~] = fileparts(pr2_file_name);
        xml_file_full_path = fullfile(xml_file_directory, [xml_file_basename '.xml']);
    end

    % Remove results folder from previous execution if still there
    results_dir = fullfile(xml_file_directory, 'results/');
    if exist(results_dir, 'dir')
        overwrite_results = questdlg('Flow results exist in working directory. Overwrite?','Warning', 'Yes','No','No');
        if strcmp(overwrite_results, 'Yes')
            command = ['rm -rf ' results_dir];
            ret = system(command);
            if ret ~= 0
                uiwait(msgbox('Cannot overwrite results folder. Aborting.', 'Flow simulation failed'));
                return
            end
        else
            return
        end
    end

    prompt = {'Number of cores to be used by HemeLB:'};
    dlg_title = 'User input';
    max_num_cores = num2str(feature('numcores'));
    default_ans = {max_num_cores};
    dialog_output = inputdlg(prompt, dlg_title, 1, default_ans);
    if isempty(dialog_output)
        % The user canceled the previous dialog box, don't continue.
        return
    end
    hemelb_num_cores = dialog_output{1};

    set(handles.wait_message,'visible','on');
    pause(0.1); % short pause ensures wait message is displayed before calling system()
    command = ['mpirun -n ' hemelb_num_cores ' hemelb -i 0 -in ' xml_file_full_path];
    status = system(command);
    set(handles.wait_message,'visible','off');

    if status==0
        title ='Flow simulation succeeded!';
        message = 'Flow simulation succeeded.';

        command = ['python PlotWSSHistogram.py ' xml_file_directory];
        system(command);
    else
        title = 'Flow simulation failed';
        message = ERROR_MESSAGE;
    end
    uiwait(msgbox(message, title));
end


%% Callbacks for third tab
%--------------------------------------------------------------------------
function [] = save_analysis_data_callback(hObject,eventdata,handles)
positionData = getappdata(handles.image_window_3,'positionData');
cellLength = getappdata(handles.image_window_3,'cellLength');
flowLength = getappdata(handles.image_window_3,'flowLength');
angles = getappdata(handles.image_window_3,'cellFlowAngle');
scalarProduct = getappdata(handles.image_window_3,'scalarProduct');
n_regions = getappdata(handles.image_window_3,'n_regions');
grp = double(getappdata(handles.image_window_3,'groups'));

% Group the data into a single matrix 
combined_data = [positionData(:,1:2),cellLength,flowLength,angles,scalarProduct,grp'];
[data_file,data_dir] = uiputfile('*.csv');
file_to_write = fullfile(data_dir,data_file);
fid = fopen(file_to_write,'w');
fprintf(fid,'# PosX,PosY,Cell Length,Flow Length,Angle,Scalar Product,Group\n');
fclose(fid);
dlmwrite(file_to_write,combined_data,'-append');

%--------------------------------------------------------------------------
function [] = show_vectors_callback(hObject,eventdata,handles)
CombineFlowAndPolarity(hObject,eventdata,handles);

%--------------------------------------------------------------------------
function [] = analyse_regions_callback(hObject,eventdata,handles)

set(handles.main_GUI,'currentaxes',handles.image_window_3);
h_mask = getappdata(handles.image_window_3,'h_mask');
% Write the polygon regions to a label mask
n_regions = getappdata(handles.image_window_3,'n_regions');
poly_region = getappdata(handles.image_window_3,'poly_region');
colour_seq = getappdata(handles.image_window_3,'colour_seq');
% setappdata(handles.image_window_3,'cmap',cmap);
% [n_rows,n_cols] = size(getappdata(handles.image_window_3,'mask_image'));
n_rows = size(getappdata(handles.image_window_3,'mask_image'),1);
n_cols = size(getappdata(handles.image_window_3,'mask_image'),2);

do_stats = true;

if n_regions < 1
    return;
end
region_matrix = zeros(n_rows,n_cols,n_regions);
% size(region)
for i = 1:n_regions
    region_matrix(:,:,i) = i*poly_region(i).createMask(h_mask);
end
setappdata(handles.tab3,'region_matrix',region_matrix);
hold off;
title('Regions')
nanregion = region_matrix;
nanregion(region_matrix==0) = nan;
label_mask = uint8(min(nanregion,[],3,'omitnan'));
label_mask(isnan(label_mask)) = 0;

overlay = ind2rgb(label_mask,colour_seq);

% Display label mask behind semi-transparent image
if isappdata(handles.image_window_3,'h_overlay_mask')
    delete(getappdata(handles.image_window_3,'h_overlay_mask'));
end
hold on;
h_overlay_mask = imshow(overlay);
% set(hMask,'AlphaData',0.25);
setappdata(handles.image_window_3,'h_overlay_mask',h_overlay_mask);

% Put the polygons back on top in case we want to move points.
for i = 1:n_regions
    uistack(h_overlay_mask,'bottom')
end
hold off;

cellFlowAngle = getappdata(handles.image_window_3,'cellFlowAngle');
flowLength = getappdata(handles.image_window_3,'flowLength');
cellLength = getappdata(handles.image_window_3,'cellLength');
scalarProduct = getappdata(handles.image_window_3,'scalarProduct');
positionData = getappdata(handles.image_window_3,'positionData');

% Assign region group according to label mask at that position
for pt = 1:size(flowLength,1)
    % Round the index to find the label designation of the nearest pixel
    grp(pt) = label_mask(round(positionData(pt,2)),round(positionData(pt,1)));
end

setappdata(handles.image_window_3,'label_mask',label_mask);
setappdata(handles.image_window_3,'n_regions',n_regions);
setappdata(handles.image_window_3,'groups',grp);

%% Run analysis on regions
% Let user decide whether to make a group from unselected region
plotBackground = 0;

if n_regions > 0
    figure(5)
    if plotBackground
        gscatter(flowLength,scalarProduct,grp,colour_seq(2:end,:));
    else
        toPlot = find(grp>0);
        gscatter(flowLength(toPlot),scalarProduct(toPlot),grp(toPlot),colour_seq(2:end,:));
    end
    title('Scatter Plot')
    flow_variable = getappdata(handles.tab3,'flow_variable');
    flow_var_str = {'WSS','Velocity','Pressure'};
    xlabel(flow_var_str{flow_variable});
    ylabel('Scalar Product');
    
    for r = 1:n_regions
        figure(10+r)
        MIN_BINS = 7;
        number_of_bins = max(MIN_BINS,round(sqrt(length(cellFlowAngle(grp==r)))));

        h_region(r) = rose_with_stats(cellFlowAngle(grp == r),number_of_bins,do_stats);
        set(h_region(r),'Color',colour_seq(r+1,:));
        title(sprintf('Angular distribution for region %d',r))
    end
end

%--------------------------------------------------------------------------
function [] = subdivide_regions_callback(hObject,eventdata,handles)
% Delete any existing regions
if isappdata(handles.image_window_3,'poly_region')
    poly_region = getappdata(handles.image_window_3,'poly_region');
    h_mask_overlay = getappdata(handles.image_window_3,'h_mask_image');
    n_regions = getappdata(handles.image_window_3,'n_regions');
    for i = 1:n_regions
        delete(poly_region(i));
        delete(h_mask_overlay);
        setappdata(handles.image_window_3,'n_regions',0);
    end
end
% Ask user how many regions
prompt = {'How many regions:'};
dlg_title = 'Enter Number of Regions';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
n_regions = str2num(answer{:});
setappdata(handles.image_window_3,'n_regions',n_regions);
SubDivideImage(handles,n_regions);

%--------------------------------------------------------------------------
function [] = get_wss_nuclei(hObject,eventdata,handles)
if getappdata(handles.tab3, 'use_last_simulation_results')
    pr2_file_name = getappdata(handles.main_GUI,'pr2_file_name');
    if isempty(pr2_file_name)
        uiwait(msgbox('No last simulation results available. Untick check box to choose manually.'));
        return;
    end
    [results_dir, ~, ~] = fileparts(pr2_file_name);
else
    results_dir = uigetdir('.', 'Choose a HemeLB simulation results folder');
    if ~exist([results_dir, '/results'],'dir')
        uiwait(msgbox('The selected folder does not contain a results/ subfolder'));
        return;
    end
end

[nucliGolgiCSVFile, nucliGolgiCSVFileFolder] = uigetfile('*.csv', 'Choose file with cell nuclei and golgi locations');
if isequal(nucliGolgiCSVFile,0)
    % User cancelled previous dialog
    return
end

dialogOutput = inputdlg('Please provide pixel/micron ratio:');
if (isempty(dialogOutput))
    % User cancelled previous dialog
    return
end
pixelsPerM = 1e6 * str2double(dialogOutput{1});
if isnan(pixelsPerM)
    uiwait(msgbox('Wrong input. Try again.'));
    return;
end

switch getappdata(handles.tab3, 'flow_variable')
    case 1
        flowVariable = 'wall_shear_stress';
    case 2
        flowVariable = 'velocity';
    case 3
        flowVariable = 'pressure';
    otherwise
        disp('Unknown flow variable when trying to call python script');
        return;
end

command = ['python ComputeWSSAtCellNuclei.py ' fullfile(nucliGolgiCSVFileFolder, nucliGolgiCSVFile) ' ' num2str(pixelsPerM) ' ' results_dir ' ' flowVariable];
ret = system(command);
if ret ~= 0
    title = ['Cannot get ' flowVariable ' at cell nuclei'];
    uiwait(msgbox(ERROR_MESSAGE, title));
    return;
end

%--------------------------------------------------------------------------
function [] = flow_variable_dropdown_callback(hObject, eventdata, handles)
setappdata(handles.tab3, 'flow_variable', hObject.Value);

%--------------------------------------------------------------------------
function [] = use_last_results_checkbox_callback(hObject, eventdata, handles)
setappdata(handles.tab3, 'use_last_simulation_results', hObject.Value);

%--------------------------------------------------------------------------
function [] = use_previous_model_setup_checkbox_callback(hObject, eventdata, handles)
setappdata(handles.tab1, 'use_previous_model_setup', hObject.Value);


%% Non-GUI functions
%--------------------------------------------------------------------------
function [] = Show_vectors(mask_image,data,handles)
figure(1)
hMask = imshow(uint8(mask_image)*55 + 200);
setappdata(handles.image_window_3,'h_mask',hMask);
set(hMask,'AlphaData',0.5);

flow_variable = getappdata(handles.tab3,'flow_variable');
flow_var_str = {'WSS','Velocity','Pressure'};
title(sprintf('%s and Cell Polarity data',flow_var_str{flow_variable}));
hold on;
quiver(data(:,1),data(:,2),data(:,3),data(:,4),0,'Color','b');
quiver(data(:,1),data(:,2),data(:,5),data(:,6),'r');
legend({'Cell polarity',flow_var_str{flow_variable}});
plot(data(:,1),data(:,2),'k.','MarkerSize',2);
axis image ij;

%--------------------------------------------------------------------------
function data = CombineFlowAndPolarity(hObject, eventdata,handles)
% Make sure we're plotting in the right place
set(handles.main_GUI,'currentaxes',handles.image_window_3);

[plexus_image_file,plexus_image_dir] = uigetfile('.tif','Choose plexus image file'); 
mask_image = im2bw(imread(fullfile(plexus_image_dir,plexus_image_file)));
setappdata(handles.image_window_3,'mask_image',mask_image);

[cell_data_file,root_dir] = uigetfile('*.csv','Choose polarity data file');
flow_data_file = strrep(cell_data_file,'.csv','_flow.csv');
xy_details = load(fullfile(root_dir,cell_data_file));

do_stats = true;

% This allows a label mask to be used to identify datapoints according to a region
% Initialised as zeros if no subdivision is required
% Note: "csvread" doesn't like blank lines, "load" doesn't seem to mind.
flowData = load(fullfile(root_dir,flow_data_file));

% Filter for very large values (typical max around 40?)
% TODO: 
%   - This should be scaled to pixel physical size
%   - Maybe pass as an argument or option via GUI?
%   - or include in a separate CONSTANTS file?
max_flow_length = inf;
flow_data_to_exclude = find(hypot(flowData(:,1),flowData(:,2)) > max_flow_length); 

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

scalarProduct = dot([data(:,3),data(:,4)],[data(:,5),data(:,6)],2);

%% BASIC PLOTS OF GLOBAL VALUES
figure(2)
num_bins_scalar_product = sqrt(length(scalarProduct));
hist(scalarProduct,num_bins_scalar_product);
title('Histogram of scalar product');

Show_vectors(mask_image,data,handles);

[cellAngle,cellLength] = cart2pol(data(:,3),data(:,4));
[flowAngle,flowLength] = cart2pol(data(:,5),data(:,6));

cellFlowAngle = cellAngle - flowAngle;  % In radians by default!
cellFlowAngleDeg = cellFlowAngle*180/pi;

figure(4)
MIN_BINS = 7;
number_of_bins = max(MIN_BINS,round(sqrt(length(cellFlowAngle))));
h_rose_all = rose_with_stats(cellFlowAngle,number_of_bins,do_stats);
title('Polar Histogram of all cell-flow angles');

% SAVE APP DATA
setappdata(handles.image_window_3,'cellFlowAngle',cellFlowAngle);
setappdata(handles.image_window_3,'flowLength',flowLength);
setappdata(handles.image_window_3,'cellLength',cellLength);
setappdata(handles.image_window_3,'scalarProduct',scalarProduct);
setappdata(handles.image_window_3,'positionData',data);

%--------------------------------------------------------------------------
function [label_mask] = SubDivideImage(handles,n_regions)
input_image = getappdata(handles.image_window_3,'mask_image');
if n_regions == 0
    label_mask = zeros(size(input_image));
    return;
end

% Use the colour sequence used in plots with black for 0.
colour_seq = [0 0 0; lines(n_regions)];
setappdata(handles.image_window_3,'colour_seq', colour_seq);

for i = 1:n_regions
    poly_region(i) = impoly;
    setappdata(handles.image_window_3,'poly_region',poly_region);
    setColor(poly_region(i),colour_seq(i+1,:));
end
