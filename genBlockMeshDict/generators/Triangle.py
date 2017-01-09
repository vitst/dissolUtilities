import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class Triangle(gbm.AbstractBaseGenerator):
  '''
  This is the generator for cylindrical meshes.
  The mesh is created with an inner rectangular prism, defined by x_rect, y_rect.
  '''

  def __init__(self):
    self.__genName__ = 'triangle'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_rect',  0.75, 'size of the rectangular prism in x direction'])
    self.parameters.append(['y_rect',  0.25, 'size of the rectangular prism in y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 20, 'number of cells in X direction'])
    self.parameters.append(['y_res', 20, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 4, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])
    self.parameters.append(['n_surf',   32, 'number of points on the edges of the triangle'])
    self.parameters.append(['beta',   1.0, 'Y direction elliptical radius'])
    self.parameters.append(['alpha',   1.5, 'X direction elliptical radius'])

  def ellipsoid_function(self, theta, delt, axis):
    alpha = self.parameters[11][1]
    beta = self.parameters[10][1]
    n_surf = self.parameters[9][1]
    #axis 0 = x, 1 = y
    if axis == 0:
    	func = alpha * np.cos(theta * np.pi/4 + delt * np.pi/2/n_surf)
    if axis == 1:
    	func = beta * np.sin(theta * np.pi/4 + delt * np.pi/2/n_surf)
    return func

  def y_ellipsoid_function(self, theta, delt):
    beta = self.parameters[10][1]
    n_surf = self.parameters[9][1]
    return beta * np.sin(theta * np.pi/4 + delt * np.pi/2/n_surf)
    

  def createBlockMeshDict(self, dictFileName):
  	# read parameters
    lines = self.read_dict(dictFileName)
    empty, x_rect = self.check_par('x_rect', lines)
    empty, y_rect = self.check_par('y_rect', lines)
    empty, Lz = self.check_par('Lz', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_G = self.check_par('x_G', lines)
    empty, y_G = self.check_par('y_G', lines)
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

    xmax = 2.0 + 3.0 * np.sqrt(3.0)
    x = np.linspace(0.0, xmax, 600)
    y = []

    for xi in x:
      yi = 0.0
      if xi<0.5:
        yi = 0.0 + np.sqrt( 1.0 - np.power( xi-1.0, 2) )
      elif xi>=0.5 and xi<= (0.5+3*np.sqrt(3.0)):
        yi = (xi - 0.5) * np.tan( np.pi / 6.0 ) + 1.0 * np.cos(np.pi/6.0)
      else:
        yi = 3.0 + np.sqrt( 1.0 - np.power( xi-1.0-3*np.sqrt(3.0), 2) )

      y.append(yi)

    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(    0,  0,  -1)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( xmax, -3,  -1)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( xmax,  3,  -1)   # 2

    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(    0,  0,  2)   # 3
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( xmax, -3,  2)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( xmax,  3,  2)   # 5

    fileCont += ');\n'

	 
    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'
    
    #******Triangle******
    fileCont += '    spline 0 1 (\n'
    for i in range(len(x)):
      fileCont += '        ({0:8g} {1:8g} {2:8g})\n'.format( x[i], -y[i],  -1)
    fileCont += '    )\n'

    fileCont += '    spline 0 2 (\n'
    for i in range(len(x)):
      fileCont += '        ({0:8g} {1:8g} {2:8g})\n'.format( x[i],  y[i],  -1)
    fileCont += '    )\n'

    fileCont += '    spline 3 4 (\n'
    for i in range(len(x)):
      fileCont += '        ({0:8g} {1:8g} {2:8g})\n'.format( x[i], -y[i],  2)
    fileCont += '    )\n'

    fileCont += '    spline 3 5 (\n'
    for i in range(len(x)):
      fileCont += '        ({0:8g} {1:8g} {2:8g})\n'.format( x[i],  y[i],  2)
    fileCont += '    )\n'
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
    fileCont += '    hex (0  1  2  0   3   4  5  3) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
    fileCont += ');\n'

    ######################################################
    # add boundaries                                     #
    ######################################################
    fileCont += '\n'
    fileCont += 'boundary\n'
    fileCont += '(\n'

    fileCont += '    aux1\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1  4  3)\n'
    fileCont += '        ( 1   2  5  4)\n'
    fileCont += '        ( 0   2  5  3)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'
    
    fileCont += '    aux2\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1  2  0)\n'
    fileCont += '        ( 3   4  5  3)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'

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
