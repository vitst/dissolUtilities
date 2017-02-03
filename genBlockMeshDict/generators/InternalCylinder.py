import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'Vitaliy Starchenko'

class InternalCylinder(gbm.AbstractBaseGenerator):
  '''
  This is the generator for cylindrical geometries, where the external wall is unsoluble and internal is soluble.
  The flow is in between two cylinders.
  The mesh is created with an inner rectangular prism, defined by x_rect, y_rect.
  '''

  def __init__(self):
    self.__genName__ = 'intcylinder'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_rect',   50.0,  'size of the rectangular prism in x direction (outer)'])
    self.parameters.append(['y_rect',   50.0,  'size of the rectangular prism in y direction (outer)'])
    self.parameters.append(['Lz',       500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['Width',    1.0,   'width of the fluid layer'])
    self.parameters.append(['tang_res', 16,    'resolution in tangential direction'])
    self.parameters.append(['norm_res', 8,     'resolution in normal direction'])
    self.parameters.append(['z_res',    512,   'resolution in Z direction'])
    self.parameters.append(['norm_G',   10.0,  'grading in normal direction'])
    self.parameters.append(['z_G',      10.0,  'grading in Z direction'])
    self.parameters.append(['n_surf',   500,    'number of points on the edges of the cylinder'])

    self.parameters.append(['beta',     1.0, 'Y direction elliptical radius'])
    self.parameters.append(['alpha',    1.5, 'X direction elliptical radius'])

  def ellipsoid_function(self, N, A, B):

    theta = np.linspace(0, np.pi, N, endpoint=True)

    x = np.zeros(N)
    y = np.zeros(N)

    for i in len(theta):
      x[i] = A * np.cos(theta[i])
      y[i] = B * np.sin(theta[i])

    return x, y

  def createBlockMeshDict(self, dictFileName):
  	 # read parameters
    lines = self.read_dict(dictFileName)
    empty, x_rect = self.check_par('x_rect', lines)
    empty, y_rect = self.check_par('y_rect', lines)
    empty, Lz = self.check_par('Lz', lines)
    empty, Width = self.check_par('Width', lines)

    empty, tang_res = self.check_par('tang_res', lines)
    empty, norm_res = self.check_par('norm_res', lines)
    empty, z_res = self.check_par('z_res', lines)

    empty, norm_G = self.check_par('norm_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    empty, n_surf = self.check_par('n_surf', lines)

    empty, beta = self.check_par('beta', lines)
    empty, alpha = self.check_par('alpha', lines)  	  
    
    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'

    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, 0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								  -y_rect, 0.0)   # 3
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, 0.0)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, 0.0)   # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), 0.0)   # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), 0.0)   # 7

    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1),  Lz)   # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1),  Lz)   # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect,  Lz)   #10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect,  Lz)   #11
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect,  Lz)   #12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect,  Lz)   #13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1),  Lz)   #14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1),  Lz)   #15

    fileCont += ');\n'

	 
    #x = alpha * np.cos(theta * phi + i * dphi)
    #y = R * np.sin(theta * phi + i * dphi)
  
  
    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'
    
    #******Cylinder******
    #Inlet
    fileCont += self.spline_function("0 1", 5, 0)
    fileCont += self.spline_function("6 7", 1, 0)
    fileCont += self.spline_function("0 6", 3, 0)
    fileCont += self.spline_function("1 7", 7, 0)
    #Outlet    
    fileCont += self.spline_function(" 8  9", 5, Lz)
    fileCont += self.spline_function("14 15", 1, Lz)
    fileCont += self.spline_function(" 8 14", 3, Lz)
    fileCont += self.spline_function(" 9 15", 7, Lz)
    #********************

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    #block 0
    fileCont += '    hex (0  1  3  2   8   9  11  10) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
    #block 1            
    fileCont += '    hex (7  6  4  5  15  14  12  13) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)

    #block 2            
    fileCont += '    hex (6  0  2  4  14   8  10  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
                
    #block 3            
    fileCont += '    hex (1  7  5  3   9  15  13  11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)       
                
    #block 3            
    fileCont += '    hex (2  3  5  4  10  11  13  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, x_res, z_res, 1, 1, z_G)          
                
    fileCont += ');\n'

    ######################################################
    # add boundaries                                     #
    ######################################################
    fileCont += '\n'
    fileCont += 'boundary\n'
    fileCont += '(\n'

    fileCont += '    outlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 8   9  11  10)\n'
    fileCont += '        (14  12  13  15)\n'
    fileCont += '        ( 8  10  12  14)\n'
    fileCont += '        ( 9  15  13  11)\n'
    fileCont += '        (10  11  13  12)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   3   2)\n'
    fileCont += '        ( 6   4   5   7)\n'
    fileCont += '        ( 0   2   4   6)\n'
    fileCont += '        ( 1   7   5   3)\n'
    fileCont += '        ( 2   3   5   4)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    walls\n'
    fileCont += '    {\n'
    fileCont += '      type wall;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   9   8)\n'
    fileCont += '        ( 6   7  15  14)\n'
    fileCont += '        ( 6   0   8  14)\n'
    fileCont += '        ( 1   7  15   9)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += ');\n'
  # *****************************************************************************
    # add boundaries                                     #
    ######################################################
    fileCont += '\n'
    fileCont += 'mergePatchPairs\n'
    fileCont += '(\n'
    fileCont += ');\n\n'

    fileCont += '// ********************************************************************* //\n'

    # writing the file
    self.writeBlockMesh( fileCont )
    return True
