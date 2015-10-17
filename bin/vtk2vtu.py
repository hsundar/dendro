#@author: Milinda Fernnado
#School of Computing University of Utah
#
#Function: Convert a series of vtk files to vtu files and write the pvtu file for render in paraview.
#

import sys
from paraview.simple import *


def vtkToVtu(file_prefix,num_proc):
    for i in range(0,num_proc):
	file_name=file_prefix+'_'+str(i)+'.vtk'
	r = LegacyVTKReader( FileNames=[file_name] )
	w = XMLUnstructuredGridWriter()
	w.FileName =file_prefix+'_'+str(i)+'.vtu'
	w.UpdatePipeline()


def pvtuFile(file_prefix,num_proc):
	filename=file_prefix+'.pvtu'
	pvtu = open(filename, 'w')
	pvtu.write('<?xml version="1.0"?> \n')
	pvtu.write('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian"> \n')
	pvtu.write('<PUnstructuredGrid GhostLevel="0"> \n')
	pvtu.write('<PPoints> \n <PDataArray type="Float32" Name="Position" NumberOfComponents="3"/> \n </PPoints> \n')
	pvtu.write('<PCells> \n   <PDataArray type="Int32" Name="connectivity" NumberOfComponents="1"/> \n')
	pvtu.write('<PDataArray type="Int32" Name="offsets"      NumberOfComponents="1"/> \n')
	pvtu.write('<PDataArray type="UInt8" Name="types"        NumberOfComponents="1"/> \n')
	pvtu.write('</PCells>\n')
	for i in range(0,num_proc):
		file_name=file_prefix+'_'+str(i)+'.vtu'	
		vtu_str='<Piece Source='+'\"'+file_name+'\"/> \n'
		pvtu.write(vtu_str)
	pvtu.write('</PUnstructuredGrid>\n')
	pvtu.write('</VTKFile>\n')

def main():
    # print command line arguments
    file_prefix=sys.argv[1]
    num_proc=int(sys.argv[2])
    ## vtk to vtu conversion
    vtkToVtu(file_prefix,num_proc)
    ## pvtu file creation
    pvtuFile(file_prefix,num_proc)
    





if __name__ == "__main__":
    main()

