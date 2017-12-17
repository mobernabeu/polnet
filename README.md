# PolNet Analysis: a software tool for the quantification of network-level endothelial cell polarity and blood flow during vascular remodelling
PolNet's graphical user interface and main pre- and post-processing algorithms. This repository contains the source code needed to build the PolNet application. You may find it easier to take a look at our [Docker containers](https://hub.docker.com/r/mobernabeu/polnet/) including a prepackaged installation of PolNet rather that building your own from source. The requirements/instructions below correspond to a PolNet installation from source.

## Requirements
* MATLAB.
* A working MATLAB [MEX](http://uk.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html) environment.
* Python including the packages: argparse, numpy, csv, matplotlib, lxml.
* A [HemeLB](https://github.com/UCL-CCS/hemelb-dev) installation. 
* A VTK installation including headers and libraries. Version 6 or higher.

## PolNet installation instructions
* Run **compile_VTP_writer.m** to compile the MEX files that provide the interface to VTK. You may need to edit the file to indicate the location of the VTK headers and libraries in your system.

## Source code being redistributed
* As part of this repository, we redistribute:
  * [CircStat2012a/](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-)
  * [voronoiSkel.m](http://uk.mathworks.com/matlabcentral/fileexchange/27543-skeletonization-using-voronoi)

## Running PolNet
* Run **PolNetGUI.m** to launch the application.
