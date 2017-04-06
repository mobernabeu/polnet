% Compile MEX files (.vtp writer and vessel crossing corrector)
switch computer
    case 'MACI64'
        % Link against a VTK 7 installation provided by Homebrew
        mex -I/usr/local/include/vtk-7.1/ -L/usr/local/lib/ -lvtkIOXML-7.1 -lvtkCommonDataModel-7.1 -lvtkCommonCore-7.1 VTP_writer.cpp
        mex -I/usr/local/include/vtk-7.1/ -L/usr/local/lib/ -lvtkIOXML-7.1 -lvtkCommonDataModel-7.1 -lvtkCommonCore-7.1 -lvtkCommonExecutionModel-7.1 -lvtkFiltersModeling-7.1 DetectAndCorrectVesselCrossings.cpp
    case 'GLNXA64'
        mex -I/usr/local/vmtk-build/Install/include/vtk-7.0/ -L/usr/local/vmtk-build/Install/lib -lvtkIOXML-7.0 -lvtkCommonDataModel-7.0 -lvtkCommonCore-7.0 VTP_writer.cpp
        mex -I/usr/local/vmtk-build/Install/include/vtk-7.0/ -L/usr/local/vmtk-build/Install/lib -lvtkIOXML-7.0 -lvtkCommonDataModel-7.0 -lvtkCommonCore-7.0 -lvtkCommonExecutionModel-7.0 -lvtkFiltersModeling-7.0 DetectAndCorrectVesselCrossings.cpp
    otherwise
        error(['I don''t know how to compile MEX files on platform ' computer])
end
