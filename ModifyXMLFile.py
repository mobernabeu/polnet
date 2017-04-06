#!/usr/bin/env python

from xml.etree import ElementTree
from argparse import ArgumentParser

def update_xml_file(filename):
    # Load automatically generated XML file
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # Get total number of timesteps
    timesteps = root.find('simulation').find('steps').attrib['value']
    
    # Add monitoring of incompressibility and convergence
    monitoring = ElementTree.SubElement(root, 'monitoring')
    ElementTree.SubElement(monitoring, 'incompressibility') 
    convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-5', 'terminate': 'false'})
    ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.001', 'units': 'm/s'})
    
    # Add definition of properties to be extracted
    extr = ElementTree.SubElement(root, 'properties')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': timesteps, 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': timesteps, 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': timesteps, 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': timesteps, 'file': 'wholegeometry-shearrate.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='shearrate')

    # Save XML file to disk
    tree.write(filename)

    
if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('xml_filename', type=str, help='XML file to be modified.')
    args = parser.parse_args()

    update_xml_file(args.xml_filename)
