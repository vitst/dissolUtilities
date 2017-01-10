import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class Triangle(gbm.AbstractBaseGenerator):
  '''
  This is the generator for triangular objects.
  The object is an isoscelles triangle, symmetric about the x axis.
  '''

  def __init__(self):
    self.__genName__ = 'triangle'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_len', 24., 'perpendicular length'])
    self.parameters.append(['y_len', 16., 'base length'])
    self.parameters.append(['x_res', 200, 'number of cells in X direction'])
    self.parameters.append(['y_res', 200, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 4, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])

  def createBlockMeshDict(self, dictFileName):
  	# read parameters
    lines = self.read_dict(dictFileName)
    empty, x_len = self.check_par('x_len', lines)
    empty, y_len = self.check_par('y_len', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_G = self.check_par('x_G', lines)
    empty, y_G = self.check_par('y_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    
    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'

    xmax = 2.0 + 3.0 * np.sqrt(3.0)
    ymax = 6.0
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

    y = np.array(y, float)
    x = x*x_len/xmax
    y = y*y_len/ymax

    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0,        0,  -1)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( x_len, -y_len/2,  -1)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( x_len,  y_len/2,  -1)   # 2

    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0,        0,  2)   # 3
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( x_len, -y_len/2,  2)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( x_len,  y_len/2,  2)   # 5

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
    fileCont += '        (0  1  4  3)\n'
    fileCont += '        (1  2  5  4)\n'
    fileCont += '        (0  2  5  3)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    aux2\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (0  1  2  0)\n'
    fileCont += '        (3  4  5  3)\n'
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
