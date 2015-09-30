import healpy as hp
import numpy as np

def readSkyMap(file):
	skymap = hp.read_map(file, verbose=False)
	npix = len(skymap)
	nside = hp.npix2nside(npix)
	sortedSky = np.sort(skymap)
	theta, phi = hp.pix2ang(nside, np.arange(0, npix)) # Construct the theta and phi arrays
	pVal = skymap[np.arange(0, npix)] # construct the pixel value array

	keep = pVal > 0. # Choose the pixels that have  non zero values [..., False, ..., True,...]
	pVal = pVal[keep] # Reset pixel value array with only non-zero elements
	order = np.argsort(-pVal) # Sort the reset pixel value array in descending order
	theta = theta[keep][order] # Sort and create the theta array values using the order of pixel value above
	phi = phi[keep][order] # Sort and create the phi value array using the order of pixel value above
	pVal = pVal[order] # Create the pixel value array using the pixel value descending order
	ra = np.rad2deg(phi) # Construct ra array
	dec = np.rad2deg(0.5*np.pi - theta) # Construct dec array
	sortedSkyPointVal = np.array([pVal, ra, dec]).transpose() # Construct the [pVal, ra, dec] array
	pValSum = np.cumsum(pVal) # Cumulative sum of pixel value array
	return ([nside, ra, dec, pVal])
