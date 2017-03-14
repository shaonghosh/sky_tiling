import numpy as np
import pylab as pl
import pickle
import sys
from math import ceil
import healpy as hp
from scipy import interpolate

import time
#from AllSkyMap_basic import AllSkyMap




class RankedTileGenerator:
	def __init__(self, skymapfile):
		self.skymap = hp.read_map(skymapfile, verbose=False)
		npix = len(self.skymap)
		self.nside = hp.npix2nside(npix)
		self.preCompDictFiles = {64:'preComputed_pixel_indices_64.dat', 
							128:'preComputed_pixel_indices_128.dat', 
							256:'preComputed_pixel_indices_256.dat',
							512:'preComputed_pixel_indices_512.dat',
							1024:'preComputed_pixel_indices_1024.dat',
							2048:'preComputed_pixel_indices_2048.dat'}
							

	def sourceTile(self, ra, dec, tiles):
		'''
		METHOD     :: This method takes the position of the injected 
					  event and returns the tile index
					  
		ra		   :: Right ascension of the source in degrees
		dec		   :: Declination angle of the source in degrees
		tiles	   :: The tile coordinate file (in the following format)
						ID	ra_center	dec_center	
						1  24.714290	-85.938460
						2  76.142860	-85.938460
						...
		'''
		tileData = np.recfromtxt(tiles, names=True)
		Dec_tile = tileData['dec_center']
		RA_tile = tileData['ra_center']
		ID = tileData['ID']
		s = np.arccos( np.sin(np.pi*dec/180.)\
			* np.sin(np.pi*Dec_tile/180.)\
			+ np.cos(np.pi*dec/180.)\
			* np.cos(np.pi*Dec_tile/180.) \
			* np.cos(np.pi*(RA_tile - ra)/180.) )
		index = np.argmin(s) ### minimum angular distance index

		return ID[index] - 1 ### Since the indexing begins with 1.


	
	def searchedArea(self, ra, dec, resolution=None):
		'''
		METHOD     :: This method takes the position of the injected 
					  event and the sky-map. It returns the searched 
					  area of the sky-map to reach to the source lo-
					  cation. The searched area constitutes both the
					  total area (sq. deg) that needed to be search-
					  ed to reach the source location, and the total
					  localization probability covered in the proce-
					  ss.
					  
		ra		   :: Right ascension of the source in degrees
		dec		   :: Declination angle of the source in degrees
		resolution :: The value of the nside, if not supplied, 
					  the default skymap is used.
		'''
		if not resolution:
			resolution = self.nside
		n = np.log(resolution)/np.log(2)
		resolution = int(2 ** round(n)) ## resolution in powers of 2
		if resolution > 2048: resolution = 2048
		if resolution < 64: resolution = 64
		filename = self.preCompDictFiles[resolution]
		File = open(filename, 'rb')
		data = pickle.load(File)
		tile_index = np.arange(len(data))
		skymapUD = hp.ud_grade(self.skymap, resolution, power=-2)
		npix = len(skymapUD)
		theta, phi = hp.pix2ang(resolution, np.arange(0, npix))
		ra_map = np.rad2deg(phi) # Construct ra array
		dec_map = np.rad2deg(0.5*np.pi - theta) # Construct dec array
		pVal = skymapUD[np.arange(0, npix)]
		order = np.argsort(-pVal)
		ra_map = ra_map[order]
		dec_map = dec_map[order]
		pVal = pVal[order]
		s = np.arccos( np.sin(np.pi*dec/180.)\
			* np.sin(np.pi*dec_map/180.)\
			+ np.cos(np.pi*dec/180.)\
			* np.cos(np.pi*dec_map/180.) \
			* np.cos(np.pi*(ra_map - ra)/180.) )
		index = np.argmin(s) ### minimum angular distance index
		coveredProb = np.sum(pVal[0:index])
		searchedArea = index*hp.nside2pixarea(resolution, degrees=True)
		return [searchedArea, coveredProb]


	
	def ZTF_RT(self, resolution=None, verbose=False):
		'''
		METHOD		:: This method returns two numpy arrays, the first
					   contains the tile indeces of ZTF and the second
					   contains the probability values of the corresp-
					   onding tiles. The tiles are sorted based on th-
					   eir probability values.
		
		resolution  :: The value of the nside, if not supplied, 
					   the default skymap is used.
		'''
		if not resolution:
			resolution = self.nside
		n = np.log(resolution)/np.log(2)
		resolution = int(2 ** round(n)) ## resolution in powers of 2
		if resolution > 2048: resolution = 2048
		if resolution < 64: resolution = 64
		if verbose: print 'Using resolution of ' + str(resolution)
		filename = self.preCompDictFiles[resolution]
		if verbose: print filename
		File = open(filename, 'rb')
		data = pickle.load(File)
		tile_index = np.arange(len(data))
		skymapUD = hp.ud_grade(self.skymap, resolution, power=-2)
		npix = len(skymapUD)
		theta, phi = hp.pix2ang(resolution, np.arange(0, npix))
		pVal = skymapUD[np.arange(0, npix)]

		allTiles_probs = []
		for ii in range(0, len(data)):
			pTile = np.sum(pVal[data[ii]])
			allTiles_probs.append(pTile)

		allTiles_probs = np.array(allTiles_probs)
		index = np.argsort(-allTiles_probs)

		allTiles_probs_sorted = allTiles_probs[index]
		tile_index_sorted = tile_index[index]

		
		return [tile_index_sorted, allTiles_probs_sorted]
		



	def plot(self, tiles, ra, dec, resolution=None, CI=0.9):
		tileData = np.recfromtxt(tiles, names=True)
		ra_center = tileData['ra_center']
		dec_center = tileData['dec_center']

		if not resolution:
			resolution = self.nside
		n = np.log(resolution)/np.log(2)
		resolution = int(2 ** round(n)) ## resolution in powers of 2
		if resolution > 2048: resolution = 2048
		if resolution < 64: resolution = 64
		print 'Using resolution of ' + str(resolution)
		skymapUD = hp.ud_grade(self.skymap, resolution, power=-2)
		hp.mollview(skymapUD)
		hp.visufunc.projplot(dec_center, ra_center,  'c.', lonlat=True)
		hp.visufunc.projplot(ra, dec,  'r*', markersize=10, lonlat=True)
		pl.show()



	def rankGalaxies2D(self, catalog, resolution=None):
		'''
		METHOD  :: This method takes as input a galaxy catalog pickle file
				   that is generated by running the createCatalog.py script.
				   The output is the IDs of the galaxies from the catalog 
				   ranked based on their localization probability.
		
		catalog	:: A pickle file which stores a 7 col numpy array with. The 
				   columns of this array are defined below:
				   		col1 : galaxy ID
				   		col2 : distance to the galaxy
				   		col3 : Declination angle of the galaxy
				   		col4 : Right ascencion of the galaxy
				   		col5 : Closest BAYESTAR Healpix pixel to the galaxy
				   		col6 : Declination angle of the closest pixel
				   		col7 : Right ascencion of the closest pixel
				   		
		resolution :: Optional argument. allows you to fix the resolution of 
					  the skymap. Currently the catalog file has only been 
					  generated for resolution of 512. Use this value.
				   
		'''

		if not resolution:
			resolution = self.nside
		n = np.log(resolution)/np.log(2)
		resolution = int(2 ** round(n)) ## resolution in powers of 2
		if resolution > 2048: resolution = 2048
		if resolution < 64: resolution = 64
		filename = self.preCompDictFiles[resolution]
		File = open(filename, 'rb')
		data = pickle.load(File)
		tile_index = np.arange(len(data))
		skymapUD = hp.ud_grade(self.skymap, resolution, power=-2)
		npix = len(skymapUD)
		theta, phi = hp.pix2ang(resolution, np.arange(0, npix))
		ra_map = np.rad2deg(phi) # Construct ra array
		dec_map = np.rad2deg(0.5*np.pi - theta) # Construct dec array
		pVal = skymapUD[np.arange(0, npix)]
		
		catalogFile = open(catalog, 'rb')
		catalogData = pickle.load(catalogFile)
		
		indices = catalogData[:,4].astype('int') ### Indices of pixels for all galaxies
		galaxy_probs = pVal[indices] ### Probability values of the galaxies in catalog
		order = np.argsort(-galaxy_probs) ### Sorting in descending order of probability
		galaxy_indices = catalogData[:,0].astype('int') ### Indices of galaxies
		ranked_galaxies = galaxy_indices[order]
		galaxy_probs = galaxy_probs[order]
		
		return [ranked_galaxies, galaxy_probs]
		
####################END OF CLASS METHODS########################
def gaussian_distribution_function(x, mu, sigma):
		'''
		METHOD	:: Creates the gaussian function corresponding to the
				   mean and standard deviation of the limiting magnitude
				   corresponding to the given time
		'''
		return 1/(np.sqrt(2*np.pi*sigma**2))*np.exp(-(x-mu)**2/(2*sigma**2))

def apparent_from_absolute_mag(absolute_mag, source_dist_parsec):
		'''
		METHOD	:: This method computes the aparent magnitude from the absolute
		magnitude model and the distance to the source.
		'''
		return absolute_mag + 5*np.log10(source_dist_parsec/10.)

def detectability(rank, time_per_tile, total_observation_time, absolute_mag, source_dist_parsec, time_data, limmag_data, error_data = None, verbose=False):
	'''
	METHOD :: This method takes as input the time allotted per tile, 
	total observation time allotted for an event, the absolute 
	magnitude model of the source and the dependence of 
	limiting magnitude on integration time, and returns either a 
	boolean numpy object[True/False] whether the source can be detected, or
	a probability of detection if the error data (sigma) for 
	the limiting magnitude is provided.
		  
	total_observation_time 	 :: The total observation time for the event
	absolute_magnitude 		 :: The absolute magnitude of the source that 
					is to be set by the model.
	source_dist_parsec 		 :: distance to source in parsecs time_data, 
					limmag_data, error_data :: The data which 
					needs to be interpolated to give limiting 
					magnitude as a function of time. error_data, 
					if provided, will allow for a detection 
					probability to be generated as output. If 
					detectability(rank, time_per_tile, total_observation_time, absolute_mag, source_dist_parsec, time_data, limmag_data, error_data = None, verbose=False)error_data is not provided, a Boolean 
					(True/False) for detection will be output. 
	'''
	### Convert to numpy object if scalar supplied
	if isinstance(time_per_tile, (np.ndarray,)) == False:
		time_per_tile = np.array(time_per_tile)
	### Check if source tile rank has been reached. Return non-detection if not
	rank_reached		= (total_observation_time/time_per_tile).astype(int)
	rank_reached_mask	= rank_reached > rank	# True means rank is reachable
	if np.all(~rank_reached_mask):	# if rank cannot be reached for 
									# any integration time	
		if verbose: print "Tile not reached in ANY allotted observation time"
		if error_data is not None:
			return np.zeros(len(time_per_tile))
		else:
			return rank_reached_mask

	### Determine limiting magnitude as a function of time, via interpolation of data
	s = interpolate.UnivariateSpline(np.log(time_data), limmag_data, k=5)
	limmag = s(np.log(time_per_tile))
	apparent_mag = apparent_from_absolute_mag(absolute_mag, source_dist_parsec)
	### If error_data is not supplied, return Boolean True/False for detection
	if error_data is None:
		depthReached = (limmag > apparent_mag)
		if np.any(depthReached) is False:
			if verbose: print "Depth not reached in ANY allotted integration time"
		return np.logical_and(depthReached, rank_reached_mask)
		### Both Depth criteria and rank criteria should be satisfied
	### If error_data is supplied, return detection probability
	else:
		s_err = interpolate.UnivariateSpline(np.log(time_data), error_data, k=5)
		mu = limmag
		sigma = s_err(np.log(time_per_tile))
		very_large_number = 1000 #proxy for +infinity
		samples = 10**5
		x = np.linspace(apparent_mag, very_large_number, samples, endpoint = True)
		### If floats are passed to the function, return the answer rightaway
		if isinstance(mu, (float, np.float, np.float64,)) == True and  isinstance(sigma, (float, np.float, np.float64,)) == True:
			y = gaussian_distribution_function(x, mu, sigma)
			return np.trapz(y,x)
		# If an array of time_per_tile had been passed
		# Note that mu and sigma are equal length arrays, as defined above
		else:
			result = []
			for ii in range(len(mu)):
				y = gaussian_distribution_function(x,mu[ii],sigma[ii])
				result.append(np.trapz(y,x))
			return np.array(result)
