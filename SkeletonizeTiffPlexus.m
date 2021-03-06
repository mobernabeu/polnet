function [skeleton,outputFileName] = SkeletonizeTiffPlexus (pixelsPerUm, tifFileName, inputImg)
    % From Martin's e-mail:
    % I've also attached the matlab script (I'll get around to reading the git info and do it properly later). It uses the matlab image processing toolbox and a script from the file-exchange called voronoiSkel (http://www.mathworks.co.uk/matlabcentral/fileexchange/27543-skeletonization-using-voronoi/content/voronoiSkel.m), which needs qhull to be installed too. It reads a tiff file (I also have a version that can pull the info straight out of Claudio's photoshop files if that's useful).
    % If the segmentation needs tweaking, there is a 'trim' parameter in voronoiSkel that can help removing spurious branches.

    % Uses voronoiSkel from fileexchange http://www.mathworks.co.uk/matlabcentral/fileexchange/27543-skeletonization-using-voronoi/content/voronoiSkel.m
    %%
    % Following line needed for running locally on R2015a (not via Docker?)
%     setenv('PATH', [getenv('PATH') ':/usr/local/bin'])
    testPlot = 0;
    if nargin < 1
        error('Not enough arguments provided');
    elseif nargin < 2
        [fileName,dirName] = uigetfile('*.tif');
        tifFileName = fullfile(dirName,fileName);
        inputImg = imread(tifFileName);
    elseif nargin < 3
        inputImg = imread(tifFileName);
    end
    [tifDir,tifName,tifExt] = fileparts(tifFileName);
    fileName = strcat(tifName,tifExt);
    if isempty(tifDir)
        dirName = pwd;
    else
        dirName = tifDir;
    end


    % Check corners of image to see if we need to invert it
    % Assumes no boundary crossing at corners!
    corner_pixels(1) = inputImg(1,1);
    corner_pixels(2) = inputImg(1,end);
    corner_pixels(3) = inputImg(end,1);
    corner_pixels(4) = inputImg(end,end);
    
    
    if sum(corner_pixels > 0) >= 3
        % Don't invert
    elseif sum(corner_pixels > 0) <= 1
        % Invert
        inputImg = imcomplement(inputImg);
    end
    
    % Turn RGB to BW, cheat and take R channel and invert it (for Black signal
    % on white background in tiff image)
    structuring_element_radius = round(1.0 * pixelsPerUm); % approx. 1 micron in pixels, could be made user-configurable
    plexusImg = imclose(~inputImg,strel('disk', structuring_element_radius));

    %%
    fprintf('Getting skeleton...');
    [voronoiSkeleton,vertices,edges] = voronoiSkel(plexusImg);
    fprintf('finished\n');
    negPlexus = ~plexusImg; % Unused? Or used in .mat file later?
    B = bwboundaries(plexusImg);
    plexusBoundary = vertcat(B{:});

    %%
    % Produce matlab skeleton by postprocessing voronoi skeleton
    skeleton = bwmorph(voronoiSkeleton,'skel','Inf');
    branchpoints = bwmorph(voronoiSkeleton,'branchpoints');

    %%
    fprintf('Preparing plexus image...');
    % [plexusR,plexusC] = find(negPlexus);
    fprintf('finished\n');
    fprintf('Getting distances...');
    % [D,I] = pdist2([plexusC plexusR],v,'euclidean','Smallest',1);
    [radiusIndex,radius] = knnsearch(plexusBoundary, vertices);
    fprintf('finished\n');

    [path,name,ext] = fileparts(fileName);
    outputFileName = fullfile(dirName,[name '.mat']);

    save(outputFileName);

    if testPlot
        %% Plot plexus with circles at voronoi vertices (proportional to vessel diameter)
        radius(radius==0) = [];
        for i = 1:length(radius)
            plot(vertices(i,2),vertices(i,1),'o','MarkerSize',radius(i)/10);
            hold on;
        end
        axis ij image;
        hold off;
    end
end
