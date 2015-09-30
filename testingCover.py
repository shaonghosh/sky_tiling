import findTiles
import readFits
from readSkyMap import readSkyMap
import numpy as np

[nside, ra,dec, pVal] = readSkyMap('/Users/shaon/Downloads/first2Years/first2yearsData/test/10116/bayestar.fits.gz')
tileBound = np.loadtxt('1x/pixelBoundary1x.dat')
object = findTiles.cover(pVal, ra, dec, tileBound, 0.95)
coverTiles = object.contourCovering()
print len(coverTiles)
