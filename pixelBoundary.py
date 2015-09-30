import numpy as np
from math import ceil
import pylab
import sys

file = open("pixelBoundary1sqDeg.dat", "w")

data = np.loadtxt('oneSquareDeg/BlackGEMtiles1sqd.asc')

BG = 1.0
count = 0
for item in data:
    pdecU = item[2] + 0.5 * np.sqrt(BG)
    pdecD = item[2] - 0.5 * np.sqrt(BG)
    praL = item[1] - 0.5 * np.sqrt(BG)/np.cos( abs(item[2])*np.pi/180. )
    praR = item[1] + 0.5 * np.sqrt(BG)/np.cos( abs(item[2])*np.pi/180. )
    if praR - praL > 360.0:
        praL = 0.0
        praR = 360.0
    pylab.plot([praL, praR], [pdecD, pdecD], color='b', linestyle='-', linewidth=2)
    pylab.hold(1)
    pylab.plot([praR, praR], [pdecD, pdecU], color='b', linestyle='-', linewidth=2)
    pylab.plot([praR, praL], [pdecU, pdecU], color='b', linestyle='-', linewidth=2)
    pylab.plot([praL, praL], [pdecU, pdecD], color='b', linestyle='-', linewidth=2)
    file.writelines( str(praL) + '\t' + str(praR) + '\t' + str(pdecD) + '\t' + str(pdecU) + '\n' )
    count += 1
    sys.stdout.write('\r')
    sys.stdout.write(str(ceil(float(count)*100/float(len(data))*100)/100. ) + ' percent done...')
    sys.stdout.flush()

file.close

pylab.xlabel('RA')
pylab.ylabel('Dec')
pylab.title('Sky tiling for BlackGEM')
pylab.show()

