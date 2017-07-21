"""
# Copyright (C) 2017 Shaon Ghosh & Deep Chatterjee
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""

from __future__ import division

import numpy as np
import healpy as hp
import time
import emcee
import pylab as pl

from utilities import tileCover
from utilities import AllSkyMap_basic


    
def getRandomPos(ra, dec, nums=1):
	'''
	Return one or more random ra and dec from the list of ras and decs supplied
	pattern of return array [RA, RA, RA, Dec, Dec, Dec] for nums=3
	'''
	return np.hstack([np.random.choice(ra, nums), np.random.choice(dec, nums)])


class PlaceTile:
	def __init__(self, skymapData, FOV, numtiles=2):
		self.skymapData = skymapData
		self.tileObj = tileCover.GetTiles(skymapData, FOV)
		self.ra_map = self.skymapData[0]
		self.dec_map = self.skymapData[1]
		self.pVal = self.skymapData[2]
		self.numtiles = numtiles
					

	def lnprior(self, skyposition): ### Basic uniform prior across the sky
		reshaped_skyposition = np.reshape(skyposition, (2, self.numtiles))
		ra_centers = reshaped_skyposition[0]
		dec_centers= reshaped_skyposition[1]
		ignore = (np.abs(ra_centers) > 360.) + (np.abs(dec_centers) > 90.)
		if np.sum(ignore):
			return -np.inf
		else:
			return 0
 

	def lnlikelihood(self, skyposition_cent):
		reshaped_skyposition = np.reshape(skyposition_cent, (2, self.numtiles))
		ra_centers = reshaped_skyposition[0]
		dec_centers= reshaped_skyposition[1]
		masked_points=np.array([]).reshape(0, 2)
		probabilitySum = 0.0
		for ii in range(self.numtiles):
			[point_ra, point_dec, p] = self.tileObj.tileCover([ra_centers[ii], dec_centers[ii]], masked=masked_points)
			points_to_be_masked = np.vstack((point_ra, point_dec)).T
			masked_points = np.vstack((masked_points, points_to_be_masked))
			probabilitySum += np.sum(p)
	
		return 4*np.log(probabilitySum)


	def lnpost(self, skypos):
		lp = self.lnprior(skypos)
		if not np.isfinite(lp):
			return -np.inf
		return lp + self.lnlikelihood(skypos)
		
		
	
	def getSamples(self):
		ndim, nwalkers = 2*self.numtiles, 4*self.numtiles
		include = np.cumsum(self.pVal) < 0.25
		include[np.sum(include)] = True
		ra_included = self.ra_map[include]
		dec_included = self.dec_map[include]
		p0 = []
		[p0.append(getRandomPos(self.ra_map, self.dec_map, nums=self.numtiles)) for _ in range(nwalkers)]
				    
		sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnpost)
		start = time.time()
		pos, prob, state = sampler.run_mcmc(p0, 100)
		print 'Burned in...'
		print 'Time taken to burn in = ' + str(time.time() - start)
		sampler.reset()
		start = time.time()
		result = sampler.run_mcmc(pos, 10000)
		end = time.time()
		print 'Acceptance fraction = ' + str(np.mean(sampler.acceptance_fraction[:]))
		print 'Time taken to finish the MCMC = ' + str(end - start)
		tileCenters = sampler.flatchain
		
		# save data in pickle file
		# if file is loaded as dat, dat[:,0] are RAs
		# dat[:,1] are Decs for the first tile and so on
		import cPickle as pickle
		print "\n\nSaving data to pickle file"
		outfile = open('samples%dTile.pickle'%(self.numtiles,), 'wb')
		data    = tileCenters.T
		output  = []
		# make the pickle file in format suggested by Shaon
		for ii in range(self.numtiles):
			output.append(data[ii])
			output.append(data[self.numtiles + ii])
		output = np.vstack(output).T
		pickle.dump(output, outfile)
		outfile.close()
            
		return output
		
		
	def optimizeBins(self, ra, dec, masked_points=np.array([]).reshape(0, 2)):
		binTrials = int(np.sqrt(len(ra)))
		bins_array = np.arange(2, binTrials + 1)
		probs_allChosenTiles = np.array([])
	
		for bins in bins_array:
			[hist, ra_bin, dec_bin] = np.histogram2d(ra, dec, bins)
			ra_bin_cent = 0.5*(ra_bin[1:] + ra_bin[:-1])
			dec_bin_cent = 0.5*(dec_bin[1:] + dec_bin[:-1])
			max_pos =np.argmax(hist)
			x = int(np.argmax(hist)/bins)
			y = np.argmax(hist)%bins
			ra_peak = ra_bin_cent[x]
			dec_peak = dec_bin_cent[y]
			[_, _, probs_chosenTile] = self.tileObj.tileCover([ra_peak, dec_peak], masked_points)
			prob_encl = np.sum(probs_chosenTile)
			probs_allChosenTiles = np.append(probs_allChosenTiles, prob_encl)
		
		bin_max = bins_array[np.argmax(probs_allChosenTiles)]

		[hist_max, ra_bin_max, dec_bin_max] = np.histogram2d(ra, dec, bin_max)
		ra_bin_cent = 0.5*(ra_bin_max[1:] + ra_bin_max[:-1])
		dec_bin_cent = 0.5*(dec_bin_max[1:] + dec_bin_max[:-1])

	
		y = np.argmax(hist_max)%bin_max
		x = int(np.argmax(hist_max)/bin_max)
		ra_peak = ra_bin_cent[x]
		dec_peak = dec_bin_cent[y]

		[ra_chosen_tile, dec_chosen_tile, probs_chosenTile] = self.tileObj.tileCover([ra_peak, dec_peak])
		masked_points = np.vstack((ra_chosen_tile, dec_chosen_tile)).T
	
		return [ra_peak, dec_peak, masked_points]
		

	def localizeTC(self, reference=None, samples=None, verbose=False):
	
		if samples is None:
			samples = self.getSamples()
		else:
			print '\n\nReading data from pickled files...'

		ra_samples = []
		dec_samples = []
		for ii in range(0, self.numtiles):
			ra_samples.append(samples[:,ii*2])
			dec_samples.append(samples[:,(ii*2)+1])
			


		masked_points=np.array([]).reshape(0, 2)
		probabilitySum = 0.0
		RA_Peak_list = []
		Dec_Peak_list = []
		for ii in range(self.numtiles):
			[ra_peak, dec_peak, points_to_be_masked] = self.optimizeBins(ra_samples[ii], dec_samples[ii], masked_points)
			masked_points = np.vstack((masked_points, points_to_be_masked))
			RA_Peak_list.append(ra_peak)
			Dec_Peak_list.append(dec_peak)

		return [RA_Peak_list, Dec_Peak_list]




class PlotterClass(PlaceTile):
	def __init__(self, ra_peaks, dec_peaks, FOV):
		self.ra_peaks = ra_peaks
		self.dec_peaks = dec_peaks
		self.FOV = FOV
		

	def plotter(self, ra_90, dec_90, filename):

		pl.figure(figsize=(80,70))
		pl.rcParams.update({'font.size': 60})
		m = AllSkyMap(projection='hammer')
		RAP_map, DecP_map = m(ra_90, dec_90) ### 90% sky localization region
		m.drawparallels(np.arange(-90.,120.,20.), color='grey', 
						labels=[False,True,True,False], labelstyle='+/-')
		m.drawmeridians(np.arange(0.,420.,30.), color='grey')
		m.drawmapboundary(fill_color='white')
		lons = np.arange(-150,151,30)
		m.label_meridians(lons, fontsize=60, vnudge=1, halign='left', hnudge=-1) 
		m.plot(RAP_map, DecP_map, 'k.', markersize=10, alpha=0.3) 

		
		for ii in range(0, len(self.ra_peaks)):
			[dec_down, dec_up,
			ra_down_left, ra_down_right, 
			ra_up_left, ra_up_right] = tileCover.getTile(self.FOV, self.ra_peaks[ii], self.dec_peaks[ii])
			
			RAP_peak, DecP_peak = m(self.ra_peaks[ii], self.dec_peaks[ii])
	
			RAP1, DecP1 = m(ra_up_left, dec_up)
			RAP2, DecP2 = m(ra_up_right, dec_up)
			RAP3, DecP3 = m(ra_down_left, dec_down)
			RAP4, DecP4 = m(ra_down_right, dec_down)

			m.plot(RAP_peak, DecP_peak, 'r.', markersize=20, mew=1)
		
			m.plot([RAP1, RAP2], [DecP1, DecP2],'r-', linewidth=4) 
			m.plot([RAP2, RAP4], [DecP2, DecP4],'r-', linewidth=4) 
			m.plot([RAP4, RAP3], [DecP4, DecP3],'r-', linewidth=4) 
			m.plot([RAP3, RAP1], [DecP3, DecP1],'r-', linewidth=4) 


		pl.savefig(filename)
# 		pl.show()

	
	









