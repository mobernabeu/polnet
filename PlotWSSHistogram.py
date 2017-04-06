#!/usr/bin/env python

import argparse
import os.path
import numpy as np
from hemeTools.parsers.extraction import ExtractedProperty
import matplotlib.pyplot as plt

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('results_folder', type=str, help='Folder containing simulation results.')
    args = parser.parse_args()
    resultsFolder =  args.results_folder

    # Reading in flow simulation results
    fileToCompare = 'surface-tractions.xtr'
    fieldToCompare = 'tangentialprojectiontraction'
    resultsFileName = os.path.join(resultsFolder, 'results', 'Extracted', 
                                   fileToCompare)
    resultsProperties = ExtractedProperty(resultsFileName)
    lastTimestep = resultsProperties.times[-1]
    resultsLastTimeStep = resultsProperties.GetByTimeStep(lastTimestep)

    wssMagnitude = np.linalg.norm(getattr(resultsLastTimeStep, fieldToCompare), axis=1)

    # Set the bin width 
    binWidth = 0.25
    bins = np.arange(min(wssMagnitude), max(wssMagnitude) + binWidth, binWidth)

    # Compute weights such that histogram appears normalised
    weights = np.ones_like(wssMagnitude)/float(len(wssMagnitude))
    
    # Plot WSS histogram
    plt.hist(wssMagnitude, bins=bins, weights=weights)
    plt.xlabel('Wall shear stress (Pa)')
    plt.ylabel('Probability')
    plt.title('Computed wall shear stress histogram')
    xMax = 15
    plt.xlim([0, xMax])
    plt.xticks(np.arange(0, xMax+1, 1))
    plt.show()

