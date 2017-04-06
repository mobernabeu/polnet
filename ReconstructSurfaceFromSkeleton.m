function [delta_x, delta_t, time_steps, final_STL_file_name] = ReconstructSurfaceFromSkeleton(filename, crossing_list, diameterInLatticeUnits, viscosityPhysicalUnits, tau, useCenterlineModeller)

    function [absolute_path] = ConstructFileAbsolutePath(dataset_path, dataset_name, dataset_suffix, file_extension)
        if strcmp(dataset_suffix, '')
            absolute_path = fullfile(dataset_path, [dataset_name '.' file_extension]);
        else
            absolute_path = fullfile(dataset_path, [dataset_name '_' dataset_suffix '.' file_extension]);
        end
    end

    if nargin < 2
        crossing_list = [];
    end

    if nargin < 3
        diameterInLatticeUnits = 7;
    end

    if nargin < 4
        viscosityPhysicalUnits = 3.85e-6;
    end

    if nargin < 5
        tau = 0.6;
    end

    if nargin <6
        useCenterlineModeller = false;
    end

    % Depending on whether we assume that the vessels remain circular in cross
    % section (1.0) or they collapse after fixation (2/pi).
    radiiFudgeFactor = 1.0;
    %radiiFudgeFactor = 2 / pi;

    % Extract datasetName to be used as identifier
    [datasetPath,datasetName,ext] = fileparts(filename);
    assert(strcmp(ext, '.mat'), 'Wrong file extension, it should be ''.mat''. Use SkeletonizeTiffPlexus to process a ''.tif'' before calling the current function.')
    
    % Load dataset
    load(filename, 'radius', 'vertices', 'edges', 'pixelsPerUm');

    pixelToUm = 1.0 / pixelsPerUm;

    % Choose according to the variable containing the radii information (see
    % comments above)
    %radiiVariable = radii;
    radiiVariable = radius;

    % Run writer
    lengthRadiusForEachEdge = VTP_writer(vertices, edges, radiiVariable, ConstructFileAbsolutePath(datasetPath, datasetName, '', 'vtp'), pixelToUm, radiiFudgeFactor);

    % Histogram of length covered by segments of a given diameter
    largestRadius = ceil(max(lengthRadiusForEachEdge(:,2)));
    range = 0+eps:0.05:largestRadius;
    lengths(1:length(range)) = 0;
    for edgeId=1:size(lengthRadiusForEachEdge,1),
        for i=1:length(range)-2,
            if (lengthRadiusForEachEdge(edgeId,2) >= range(i) && lengthRadiusForEachEdge(edgeId,2) < range(i+1)),
                lengths(i) = lengths(i) + lengthRadiusForEachEdge(edgeId,1);
            end
        end
        i = length(range)-1;
        if (lengthRadiusForEachEdge(edgeId,2) >= range(i)),
            lengths(i) = lengths(i) + lengthRadiusForEachEdge(edgeId,1);
        end
    end
    
    % Do some stats
    parmhat = lognfit(2*range, [], [], lengths);    
    mu = parmhat(1);
    sigma = parmhat(2);
    mean = exp(mu + sigma^2/2);
    variance = exp(2*mu + sigma^2)*(exp(sigma^2) - 1);
    mode = exp(mu - sigma^2);
    median = exp(mu);

    % Plot histogram
    figure; bar(2*range,lengths)
    xlim([0 inf])
    xlabel('Diameter (microns)')
    ylabel('Total length covered (microns)')
    title('Network diameter histogram')

    % Plot a vertical line to mark the mean of the lognormal distribution
    hold on 
    line = plot([mean mean], ylim, 'r');
    legend_txt= sprintf('mean = %.2f microns', mean);
    legend(line, legend_txt);

    % Force histogram to show before going ahead
    drawnow
        
    % Compute the value of delta x required to satisfy that no more
    % than fractionToBeCovered of the total length of the network has fewer
    % than diameterInLatticeUnits lattice sites accross.
    totalLength = sum(lengths);
    partialSum = 0;
    fractionToBeCovered = 0.05;
    for i=1:length(lengths),
        partialSum = partialSum + lengths(i);
        if partialSum > (fractionToBeCovered*totalLength),
            sprintf('%f length covered at %f um', fractionToBeCovered, range(i))
            delta_x = 2*range(i)/(diameterInLatticeUnits-1);
            sprintf('Set delta x to %f um', delta_x)
            break        
        end
    end

    % based on the fact that we want viscosityPhysicalUnits to map to the
    % value of tau specified by the user (or its default).
    delta_t = (delta_x*1e-6)^2 * (tau - 0.5) / (3 * viscosityPhysicalUnits);
    sprintf('Set delta t to %e s', delta_t)

    t_diff = (2 * largestRadius * 1e-6)^2 / viscosityPhysicalUnits;
    time_steps = t_diff/delta_t;
    sprintf('Run %f timesteps', time_steps)

    DetectAndCorrectVesselCrossings(ConstructFileAbsolutePath(datasetPath, datasetName, '', 'vtp'), ConstructFileAbsolutePath(datasetPath, datasetName, 'corrected', 'vtp'), crossing_list)
    datasetName =  [datasetName '_corrected'];

    % Using vmtkcenterlinemodeller to generate the surface from the skeleton takes forever. I checked with Luca 
    % and vmtkpolyballmodeller is a valid option as long as the network nodes are quite close together (which is our case)
    ifile = ConstructFileAbsolutePath(datasetPath, datasetName, '', 'vtp');
    ofile = ConstructFileAbsolutePath(datasetPath, datasetName, 'tubed', 'stl');
    if useCenterlineModeller
        vmtkPipeline = ['vmtk vmtkcenterlinemodeller -ifile ' ifile ' -radiusarray Radius -dimensions 256 256 64 --pipe vmtkmarchingcubes -ofile ' ofile]
    else
        vmtkPipeline = ['vmtk vmtkpolyballmodeller -ifile ' ifile ' -radiusarray Radius -dimensions 1024 1024 256 --pipe vmtkmarchingcubes -ofile ' ofile]
    end
    ret = system(vmtkPipeline);

    msg_title = 'Surface reconstruction failed';
    message = 'Please check the PolNet terminal for error messages and refer to the Troubleshooting section of the paper.';
    if ret ~= 0
        uiwait(msgbox(message, msg_title));
        final_STL_file_name = '';
        return;
    end

    % Filter configuration values 0.1 and 30 taken from VMTK tutorial "PREPARE A SURFACE FOR MESH GENERATION"
    ifile = ConstructFileAbsolutePath(datasetPath, datasetName, 'tubed', 'stl');
    final_STL_file_name = ConstructFileAbsolutePath(datasetPath, datasetName, 'tubed_smoothed', 'stl');
    vmtkPipeline = ['vmtk vmtksurfacesmoothing -ifile ' ifile ' -passband 0.1 -iterations 30 -ofile ' final_STL_file_name]
    ret = system(vmtkPipeline);

    if ret ~= 0
        uiwait(msgbox(message, msg_title));
        return;
    end
end
