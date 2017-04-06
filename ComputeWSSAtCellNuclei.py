#!/usr/bin/env python

import argparse
import os.path
import numpy as np
import csv
from hemeTools.parsers.extraction import ExtractedProperty
from pylab import plot, savefig

class BruteForceMapper(object):
    def __init__(self, nucleiLocation, searchTolerance, resultsGridPosition):
        self.nucleiLocation = nucleiLocation
        self.searchTolerance = searchTolerance
        self.resultsGridPosition = resultsGridPosition

    def MapLocationFromFineToCoarseGrid(self):
        fineGridIndices = []
        subsetMapped = []

        for iter, nucleusLocation in enumerate(self.nucleiLocation):
            if (iter % 100 == 0): 
                print "Classified {} of {}".format(iter, 
                                                   len(self.nucleiLocation))

            minIndex = np.argmin(np.linalg.norm(self.resultsGridPosition - nucleusLocation, axis=1))
            minDistance = np.linalg.norm(self.resultsGridPosition[minIndex] - nucleusLocation)

            if minDistance < self.searchTolerance:
                fineGridIndices.append(minIndex)
                subsetMapped.append(iter)

        print '{} locations could not be mapped'.format(
            len(self.nucleiLocation) - len(subsetMapped))

        return (fineGridIndices, subsetMapped)
        

    def MapSolutionFromFineToCoarseGrid(self, resultsSolution):
        assert len(resultsSolution) == len(self.resultsGridPosition)

        fineGridIndices, subsetMapped = self.MapLocationFromFineToCoarseGrid()
        plot(self.resultsGridPosition[fineGridIndices, 0], self.resultsGridPosition[fineGridIndices, 1], 'g,')

        return (resultsSolution[fineGridIndices], subsetMapped)

# Maps variable name (as inputted to the script) to a pair made of the file name
# and the name of field inside the file containing the relavant data
flowVariablesFiles = {'wall_shear_stress' : ('surface-tractions.xtr', 'tangentialprojectiontraction'),
                      'velocity' : ('wholegeometry-velocity.xtr', 'velocity'),
                      'pressure' : ('surface-pressure.xtr', 'pressure')}

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('nuclei_golgi_file', type=str, help='Name of the file with cell nuclei and golgi coordinates to be used.')
    parser.add_argument('pixel_per_m', type=str, help='Pixels per metre in original image.')
    parser.add_argument('results_folder', type=str, help='Folder containing simulation results.')
    parser.add_argument('flow_variable', type=str, choices=flowVariablesFiles.keys(), help='Flow variable to be extracted.')
    parser.add_argument('--fast', help='Uses a faster but inexact tree-based search algorithm.')

    args = parser.parse_args()

    nucleiCoordsFilename = args.nuclei_golgi_file
    datasetName = os.path.splitext(nucleiCoordsFilename)[0]
    resultsFolder =  args.results_folder
    pixelsPerMetre = float(args.pixel_per_m)
    swapXYAxes =  True
    dropGolgiCoords = True

    useBruteForceSearch = not args.fast

    # Reading in cell nuclei coordinates
    nucleiCoords = []
    with open(nucleiCoordsFilename, 'r') as csvFile:
        for line in csvFile:
            coords = line.rstrip().split(',')
            if dropGolgiCoords:
                coords = coords[:2]
            if swapXYAxes:
                # Swap X and Y coordinates to match flow model convention, if required
                coords.reverse()
            coords.append('0')
            assert(len(coords)==3)
            nucleiCoords.append([float(coord)/pixelsPerMetre for coord in coords])
    nucleiCoords = np.array(nucleiCoords)

    # Reading in flow simulation results
    fileToCompare, fieldToCompare = flowVariablesFiles[args.flow_variable]
    resultsFileName = os.path.join(resultsFolder, 'results', 'Extracted', 
                                   fileToCompare)
    resultsProperties = ExtractedProperty(resultsFileName)
    lastTimestep = resultsProperties.times[-1]
    resultsLastTimeStep = resultsProperties.GetByTimeStep(lastTimestep)

    # Plot cell nuclei location against flow model domain to check they overlap
    plot(resultsLastTimeStep.position[:,0], resultsLastTimeStep.position[:,1], 'b,')
    plot(nucleiCoords[:,0], nucleiCoords[:,1], 'r+')

    if fieldToCompare in ['tangentialprojectiontraction', 'pressure']:
        # Ignore z component of the lattices. This is because WSS and pressure
        # are obtained at the model surface and the sample points are assumed to
        # lie in the z=0 plane.
        resultsLastTimeStep.grid[:,2] = 0
        resultsLastTimeStep.position[:,2] = 0
        # Tolerance in search operation (times voxel size)
        toleranceFactor = 1
    elif fieldToCompare == 'velocity':
        # Tolerance in search operation (times voxel size)
        toleranceFactor = 1
    else:
        raise NotImplementedError('Unknown flow variable requested')

    if useBruteForceSearch:
        gridMapper = BruteForceMapper(nucleiCoords,
                                      resultsProperties.voxelSizeMetres*toleranceFactor,
                                      resultsLastTimeStep.position)
    else:
        raise NotImplementedError('Will port the non brute force approach when required')

        # Create a mapper object between the cell nuclei location and the flow model
        from PlotConvergenceAnalysis import DifferentVoxelSizeGridMapper
        gridMapper = DifferentVoxelSizeGridMapper(nucleiCoords,
                                                  resultsProperties.voxelSizeMetres*toleranceFactor, # Using this as tolerance
                                                  resultsLastTimeStep.grid,
                                                  resultsProperties.voxelSizeMetres,
                                                  resultsProperties.originMetres)

    def PassInRightNPShape(results, field):
        ret = getattr(resultsLastTimeStep, fieldToCompare)
        if len(ret.shape) == 2:
            assert ret.shape[1] == 3
            return ret
        else:
            assert len(ret.shape) == 1
            return ret[:,np.newaxis]

    # Map the flow simulation solution onto the cell nuclei location
    mappedSolution, subsetMapped = gridMapper.MapSolutionFromFineToCoarseGrid(
         PassInRightNPShape(resultsLastTimeStep, fieldToCompare))

    solutionMappedOnNucleiLocation = np.zeros((len(nucleiCoords), mappedSolution.shape[1]))
    solutionMappedOnNucleiLocation[subsetMapped] = mappedSolution

    # Write out WSS at nearby surface point to cell nuclei location
    with open(os.path.join(resultsFolder, datasetName + '_flow.csv'), 'w') as csvfile:
        csvWriter = csv.writer(csvfile)
        for wss in solutionMappedOnNucleiLocation:
            if len(wss) == 3:
                if swapXYAxes:
                    # Swap X and Y wss coordinates back to match Martin's convention
                    csvWriter.writerow([wss[1], wss[0]])
                else:
                    csvWriter.writerow(wss[0:2])
            else:
                csvWriter.writerow(wss)

    # Save plot with nuclei location and flow domain locations mapped to them
    savefig(os.path.join(resultsFolder, datasetName + '.png'))
