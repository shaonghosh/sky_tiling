import numpy as np
import healpy as hp
import argparse
import pickle
import time

from astropy.utils.console import ProgressBar

# tileFile = 'ZTF_tiles_set2_nowrap_indexed_new.dat'
# surrogate_fitsfile = '/Users/ghosh4/Downloads/2016_fits/288172/bayestar.fits.gz'
# target_nside = 64


parser = argparse.ArgumentParser()
parser.add_argument("-T", "--tileFile", action="store", help="Name of the tile center file")
parser.add_argument("-R", "--resolution", action="store", type=int, help="Resoultion of the sky map")
parser.add_argument("-M", "--map", action="store", help="Surrogate sky map to get pixel locations")
parser.add_argument("-o", "--outTag", action="store", help="Tag for the precomputed files")
args = parser.parse_args()


surrogate_fitsfile = args.map
target_nside = args.resolution
tileFile = args.tileFile

start = time.time()
skymap = hp.read_map(surrogate_fitsfile)
npix = len(skymap)
nside = hp.npix2nside(npix)

if nside != target_nside:
	print 'Downsampling to NSIDE = ' + str(target_nside)

skymap = hp.ud_grade(skymap, target_nside, power=-2) ### Downsampling the sky-map
npix = len(skymap)
nside = hp.npix2nside(npix)
sortedSky = np.sort(skymap)
theta, phi = hp.pix2ang(nside, np.arange(0, npix)) # Construct the theta and phi arrays
pVal = skymap[np.arange(0, npix)] # construct the pixel value array
ra = np.rad2deg(phi) # Construct ra array
dec = np.rad2deg(0.5*np.pi - theta) # Construct dec array
pixelIndex = np.arange(npix)

data = np.recfromtxt(tileFile, names=True)
RA_tile = data['ra_center'] ### RA value of the telescope fields
Dec_tile = data['dec_center'] ### Dec values of the telescope fields
tile_index = data['ID']-1 ### Indexing the tiles

closestTileIndex = []

with ProgressBar(len(pixelIndex)) as bar:
	for ii in pixelIndex:
		s = np.arccos( np.sin(np.pi*dec[ii]/180.)\
			* np.sin(np.pi*Dec_tile/180.)\
			+ np.cos(np.pi*dec[ii]/180.)\
			* np.cos(np.pi*Dec_tile/180.) \
			* np.cos(np.pi*(RA_tile - ra[ii])/180.) )
		index = np.argmin(s) ### minimum angular distance index
		closestTileIndex.append(tile_index[index])
		bar.update()

closestTileIndex = np.array(closestTileIndex)
uniqueTiles = np.unique(closestTileIndex)

pixelsInTile = []
for tile in uniqueTiles:
	whereThisTile = tile == closestTileIndex ### pixels indices in this tile
	pixelsInTile.append(pixelIndex[whereThisTile])

end = time.time()

print 'Total time taken = ' + str(end - start)

File = open(args.outTag + '_' + str(target_nside) + '.dat', 'wb')
pickle.dump(pixelsInTile, File)
File.close()

  


















