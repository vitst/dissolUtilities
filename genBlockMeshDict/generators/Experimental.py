import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'vstar'


class Experimental(gbm.AbstractBaseGenerator):
  '''
  This generator creates a geometry similar to experimental setup used for gypsum.
  '''
  def __init__(self):
    self.__genName__ = 'experimental'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description
    self.parameters.append(['Lx',  360.0, 'length of the lab setup (X direction)'])
    self.parameters.append(['Ly',    4.0, 'width of the lab setup (Y direction)'])
    self.parameters.append(['Lz',    1.0, 'initial aperture (Z direction)'])

    self.parameters.append(['Lgap', 10.0, 'a size of the inlet and the outlet'])

    self.parameters.append(['x_res', 720, 'total number of cells in X direction'])
    self.parameters.append(['y_res',   4, 'number of cells in Y direction'])
    self.parameters.append(['z_res',  16, 'number of cells in Z direction'])

    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])

  def createBlockMeshDict(self, dictFileName):
    # read parameters
    lines = self.read_dict(dictFileName)
    empty, Lx = self.check_par('Lx', lines)
    empty, Ly = self.check_par('Ly', lines)
    empty, Lz = self.check_par('Lz', lines)

    empty, Lgap = self.check_par('Lgap', lines)

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

    # Bottom plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(   -Lgap, -Ly/2., 0)     # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0.0, -Ly/2., 0)     # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(      Lx, -Ly/2., 0)     # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx+Lgap, -Ly/2., 0)     # 3

    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx+Lgap,  Ly/2., 0)     # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(      Lx,  Ly/2., 0)     # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0.0,  Ly/2., 0)     # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(   -Lgap,  Ly/2., 0)     # 7

    # Top plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(   -Lgap, -Ly/2., 1)     # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0.0, -Ly/2., 1)     # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(      Lx, -Ly/2., 1)     # 10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx+Lgap, -Ly/2., 1)     # 11

    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx+Lgap,  Ly/2., 1)     # 12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(      Lx,  Ly/2., 1)     # 13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(     0.0,  Ly/2., 1)     # 14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(   -Lgap,  Ly/2., 1)     # 15

    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'
    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    # by default symmetric grading across the fracture
    #y_Gs = '((0.5 0.5 {0:g}) (0.5 0.5 {1:g}))'.format(y_G, 1.0/y_G)
    y_Gs = '{0:g}'.format(y_G)

    # TODO: implement more general resolution calculation
    # TODO: taking into account grading
    x_res1 = int( 2*Lgap / (Lx+2*Lgap) * x_res )
    y_res1 = y_res
    z_res1 = z_res

    x_res2 = x_res - x_res1
    y_res2 = y_res
    z_res2 = z_res

    if x_res1 % 2 > 0:
      x_res1 -= 1
    x_res1 = int(x_res1/2)

    fileCont += '    hex (0  1  6  7  8  9  14  15) ({0:d} {1:d} {2:d})\n' \
                '    simpleGrading (\n' \
                '      {3:g}\n' \
                '      {4:s}\n' \
                '      {5:g}\n' \
                '    )\n'. \
                format(x_res1, y_res1, z_res1,   1, y_Gs, z_G)

    fileCont += '    hex (1  2  5  6  9  10  13  14) ({0:d} {1:d} {2:d})\n' \
                '    simpleGrading (\n' \
                '      {3:g}\n' \
                '      {4:s}\n' \
                '      {5:g}\n' \
                '    )\n'. \
                format(x_res2, y_res2, z_res2, x_G, y_Gs, z_G)

    fileCont += '    hex (2  3  4  5  10  11  12  13) ({0:d} {1:d} {2:d})\n' \
                '    simpleGrading (\n' \
                '      {3:g}\n' \
                '      {4:s}\n' \
                '      {5:g}\n' \
                '    )\n'. \
                format(x_res1, y_res1, z_res1,   1, y_Gs, z_G)

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
    fileCont += '        (10   11   12   13)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (8   9   14  15)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    X1\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (0   7   15   8)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'


    fileCont += '    X2\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (3   4   12   11)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'



    fileCont += '    solubleWall\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0  1  6  7)\n'
    fileCont += '        ( 1  2  5  6)\n'
    fileCont += '        ( 2  3  4  5)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    insolubleWall\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 9  10  13  14)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'


    fileCont += '    periodicx1\n'
    fileCont += '    {\n'
    fileCont += '      type cyclic;\n'
    fileCont += '      neighbourPatch periodicx2;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   9  8)\n'
    fileCont += '        ( 1   2  10  9)\n'
    fileCont += '        ( 2   3  11 10)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    periodicx2\n'
    fileCont += '    {\n'
    fileCont += '      type cyclic;\n'
    fileCont += '      neighbourPatch periodicx1;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 7    6   14   15)\n'
    fileCont += '        ( 6    5   13   14)\n'
    fileCont += '        ( 5    4   12   13)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += ');\n'

    ######################################################
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

