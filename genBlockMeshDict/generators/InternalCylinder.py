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
    self.parameters.append(['tang_res', 50,    'resolution in tangential direction'])
    self.parameters.append(['norm_res', 8,     'resolution in normal direction'])
    self.parameters.append(['z_res',    400,   'resolution in Z direction'])
    self.parameters.append(['norm_G',   0.1,  'grading in normal direction'])
    self.parameters.append(['z_G',      1.0,  'grading in Z direction'])
    self.parameters.append(['n_surf',   500,    'number of points on the edges of the cylinder'])

  def ellipsoid_function(self, N, A, B):

    theta = np.linspace(0, np.pi/2., N, endpoint=True)

    x = np.zeros(N)
    y = np.zeros(N)

    for i in range(N):
      x[i] = A * np.cos(theta[i])
      y[i] = B * np.sin(theta[i])

    return x, y
  
  def edge(self, n_surf, pids, x, y, currentZ):
    fileCont2 = '\n'
    fileCont2 += '    spline {} (\n'.format(pids[0])
    for i in range(n_surf):
      currentX = -x[i]
      currentY = -y[i]
      fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont2 += '    )\n'
  
    fileCont2 += '\n'
    fileCont2 += '    spline {} (\n'.format(pids[1])
    for i in range(n_surf):
      currentX =  x[n_surf - 1 - i]
      currentY = -y[n_surf - 1 - i]
      fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont2 += '    )\n'
  
    fileCont2 += '\n'
    fileCont2 += '    spline {} (\n'.format(pids[2])
    for i in range(n_surf):
      currentX = x[i]
      currentY = y[i]
      fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont2 += '    )\n'
  
    fileCont2 += '\n'
    fileCont2 += '    spline {} (\n'.format(pids[3])
    for i in range(n_surf):
      currentX = -x[n_surf - 1 - i]
      currentY =  y[n_surf - 1 - i]
      fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont2 += '    )\n'
  
    return fileCont2

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

    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'

    xO = x_rect/2.0 + Width / 2.0
    xI = x_rect/2.0 - Width / 2.0
    yO = y_rect/2.0 + Width / 2.0
    yI = y_rect/2.0 - Width / 2.0

    # Lowest plane
    # outer
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( -xO, 0.0, 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0, -yO, 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(  xO, 0.0, 0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0,  yO, 0.0)   # 3
    # inner
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( -xI, 0.0, 0.0)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0, -yI, 0.0)   # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(  xI, 0.0, 0.0)   # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0,  yI, 0.0)   # 7

    # Highest plane
    # outer
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( -xO, 0.0, Lz)   # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0, -yO, Lz)   # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(  xO, 0.0, Lz)   # 10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0,  yO, Lz)   # 11
    # inner
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( -xI, 0.0, Lz)   # 12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0, -yI, Lz)   # 13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(  xI, 0.0, Lz)   # 14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 0.0,  yI, Lz)   # 15

    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'
    
    #******Outer Cylinder******
    x, y = self.ellipsoid_function(n_surf, xO, yO)
    #Inlet
    pids = ['0 1', '1 2', '2 3', '3 0']
    fileCont += self.edge(n_surf, pids, x, y, 0.0)
    #Outlet
    pids = ['8 9', '9 10', '10 11', '11 8']
    fileCont += self.edge(n_surf, pids, x, y, Lz)
    #**************************
    
    #******Inner Cylinder******
    x, y = self.ellipsoid_function(n_surf, xI, yI)
    #Inlet
    pids = ['4 5', '5 6', '6 7', '7 4']
    fileCont += self.edge(n_surf, pids, x, y, 0.0)
    #Outlet
    pids = ['12 13', '13 14', '14 15', '15 12']
    fileCont += self.edge(n_surf, pids, x, y, Lz)
    #**************************

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    #block 0
    fileCont += '    hex (0  1  5  4   8   9  13  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(tang_res, norm_res, z_res, 1.0, norm_G, z_G)
    #block 1            
    fileCont += '    hex (1  2  6  5  9  10  14  13) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(tang_res, norm_res, z_res, 1.0, norm_G, z_G)

    #block 2            
    fileCont += '    hex (2  3  7  6  10   11  15  14) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(tang_res, norm_res, z_res, 1.0, norm_G, z_G)
                
    #block 3            
    fileCont += '    hex (3  0  4  7  11  8  12  15) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(tang_res, norm_res, z_res, 1.0, norm_G, z_G)
                
                
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
    fileCont += '        ( 8   9  13  12)\n'
    fileCont += '        ( 9  10  14  13)\n'
    fileCont += '        (10  11  15  14)\n'
    fileCont += '        (11   8  12  15)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   5   4)\n'
    fileCont += '        ( 1   2   6   5)\n'
    fileCont += '        ( 2   3   7   6)\n'
    fileCont += '        ( 3   0   4   7)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    insolubleSurf\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   9   8)\n'
    fileCont += '        ( 1   2  10   9)\n'
    fileCont += '        ( 2   3  11  10)\n'
    fileCont += '        ( 3   0   8  11)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    solubleSurf\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 4   5  13  12)\n'
    fileCont += '        ( 5   6  14  13)\n'
    fileCont += '        ( 6   7  15  14)\n'
    fileCont += '        ( 7   4  12  15)\n'
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
