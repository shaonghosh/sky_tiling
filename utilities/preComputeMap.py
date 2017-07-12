import numpy as np
import healpy as hp
import pickle
import time

from astropy.utils.console import ProgressBar


def preComputeMap(tileFile, telescope, target_nside=256):
	surrogate_fitsfile = 'utilities/bayestar.fits.gz'
	
	start = time.time()
	skymap = hp.read_map(surrogate_fitsfile, verbose=False)
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

	print '\nCreating pixel tile map for skymaps of nside = ' + str(target_nside)
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

	tile_map_file_name = 'tile_pixel_maps/preComputed_' + telescope + '_pixel_indices_' + str(target_nside) + '.dat'
	File = open(tile_map_file_name, 'wb')
	pickle.dump(pixelsInTile, File)
	File.close()
	return tile_map_file_name

