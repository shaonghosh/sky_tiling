import numpy as np
import pylab as pl
import healpy as hp
import pandas as pd
from astropy.time import Time
import ConfigParser
import rankedTilesGenerator
from astropy.table import Table
import associateBNSEvents
from AllSkyMap_basic import AllSkyMap

import os

def getTileBounds(FOV, ra_cent, dec_cent):
    dec_down = dec_cent - 0.5*np.sqrt(FOV)
    dec_up = dec_cent + 0.5*np.sqrt(FOV)

    ra_down_left = ra_cent - 0.5*(np.sqrt(FOV)/(np.cos(dec_down*(np.pi/180.))))
    ra_down_right = ra_cent + 0.5*(np.sqrt(FOV)/(np.cos(dec_down*(np.pi/180.))))
    ra_up_left = ra_cent - 0.5*(np.sqrt(FOV)/(np.cos(dec_up*(np.pi/180.))))
    ra_up_right = ra_cent + 0.5*(np.sqrt(FOV)/(np.cos(dec_up*(np.pi/180.))))
    
    return([dec_down, dec_up, ra_down_left, ra_down_right, ra_up_left, ra_up_right])

def readSkymap(skymapFile, CI=0.9):
	skymap = hp.read_map(skymapFile, verbose=False)
	skymap = hp.ud_grade(skymap, 256, power=-2)
	npix = len(skymap)
	nside = hp.npix2nside(npix)
	theta, phi = hp.pix2ang(nside, np.arange(0, npix))
	ra = np.rad2deg(phi)
	dec = np.rad2deg(0.5*np.pi - theta)
	pVal = skymap[np.arange(0, npix)]
	order = np.argsort(-pVal)
	ra = ra[order]
	dec = dec[order]
	pVal = pVal[order]
	include = np.cumsum(pVal) < CI
	include[np.sum(include)] = True
	ra_CI = ra[include]
	dec_CI = dec[include]
	pVal_CI = pVal[include]
	return [ra_CI, dec_CI, pVal_CI]
	
class MakeVideo:
	def __init__(self, skymapFile, df, configfile, FOV, event=None):
		self.FOV = FOV
		self.df = df
		if event:
			self.event=event
		else:
			self.event=None
		self.configParser = ConfigParser.ConfigParser()
		self.configParser.read(configfile)
		[self.ra_CI, self.dec_CI, self.pVal_CI] = readSkymap(skymapFile)

	def plotter(self, thisTileindex, previousTileIndices, title=None, tag=None):
		'''
		This function takes as input the indices of the tiles that we want to plot, 
		makes the the plot and returns tile index that it has just plotted.
		The "thisTileIndex" is plotted in red and the previous tiles are plotted 
		in black. If the event is an injection (or the location is known), then
		the tiles that has the event is plotted in blue.
		'''
		pl.clf()
		pl.figure(figsize=(60,40))
		pl.rcParams.update({'font.size': 70})
		pl.figure(figsize=(60,40))
		if title: pl.title(title)
		m = AllSkyMap(projection='hammer')
		if self.event: RAP_event, DecP_event = m(self.event[0], self.event[1])
		m.drawparallels(np.arange(-90.,120.,20.), color='grey',\
						labels=[False,True,True,False], labelstyle='+/-')
		m.drawmeridians(np.arange(0.,420.,30.), color='grey')
		m.drawmapboundary(fill_color='white')
		lons = np.arange(-150,151,30)
		m.label_meridians(lons, fontsize=70, vnudge=1, halign='left', hnudge=-1)
		RAP_map, DecP_map = m(self.ra_CI, self.dec_CI)
		m.plot(RAP_map, DecP_map, color='y', marker='.', linewidth=0, markersize=4, alpha=1.0) 
		if self.event: m.plot(RAP_event, DecP_event, color='b', marker='*', linewidth=0, markersize=20, alpha=1.0) 
		tileFile = self.configParser.get('tileFiles', 'tileFile')
		tileData = np.recfromtxt(tileFile, names=True)
		Dec_tile = tileData['dec_center']
		RA_tile = tileData['ra_center']
		ID = tileData['ID']	
		thisTileRA = RA_tile[thisTileindex]
		thisTileDec = Dec_tile[thisTileindex]
		[dec_down_this, dec_up_this,\
		ra_down_left_this, ra_down_right_this,\
		ra_up_left_this, ra_up_right_this] = getTileBounds(self.FOV, thisTileRA, thisTileDec)

		RAP1_this, DecP1_this = m(ra_up_left_this, dec_up_this)
		RAP2_this, DecP2_this = m(ra_up_right_this, dec_up_this)
		RAP3_this, DecP3_this = m(ra_down_left_this, dec_down_this)
		RAP4_this, DecP4_this = m(ra_down_right_this, dec_down_this)
		lw = 4
		m.plot([RAP1_this, RAP2_this], [DecP1_this, DecP2_this],'r-', linewidth=lw*1.5) 
		m.plot([RAP2_this, RAP4_this], [DecP2_this, DecP4_this],'r-', linewidth=lw*1.5) 
		m.plot([RAP4_this, RAP3_this], [DecP4_this, DecP3_this],'r-', linewidth=lw*1.5)
		m.plot([RAP3_this, RAP1_this], [DecP3_this, DecP1_this],'r-', linewidth=lw*1.5) 
		if len(previousTileIndices)>0:
			for index in previousTileIndices:
				previousTileRA = RA_tile[index]
				previousTileDec = Dec_tile[index]
				[dec_down_previous, dec_up_previous, ra_down_left_previous,\
				ra_down_right_previous, ra_up_left_previous,\
				ra_up_right_previous] = getTileBounds(self.FOV, previousTileRA, previousTileDec)
				RAP1_previous, DecP1_previous = m(ra_up_left_previous, dec_up_previous)
				RAP2_previous, DecP2_previous = m(ra_up_right_previous, dec_up_previous)
				RAP3_previous, DecP3_previous = m(ra_down_left_previous, dec_down_previous)
				RAP4_previous, DecP4_previous = m(ra_down_right_previous, dec_down_previous)
				m.plot([RAP1_previous, RAP2_previous], [DecP1_previous, DecP2_previous],'k-', linewidth=lw*1, alpha=0.5)
				m.plot([RAP2_previous, RAP4_previous], [DecP2_previous, DecP4_previous],'k-', linewidth=lw*1, alpha=0.5)
				m.plot([RAP4_previous, RAP3_previous], [DecP4_previous, DecP3_previous],'k-', linewidth=lw*1, alpha=0.5)
				m.plot([RAP3_previous, RAP1_previous], [DecP3_previous, DecP1_previous],'k-', linewidth=lw*1, alpha=0.5)

		previousTileIndices = np.append(previousTileIndices, thisTileindex)
		filenametag = self.configParser.get('plot', 'filenametag')
		if tag: filenametag = filenametag + '_' + tag
		pl.savefig('schedule_output/skyTiles_' + filenametag + '.png')
		return previousTileIndices




	def makeVideo(self, event=None, fr=4, tag=None, keepdata=False):
		'''
		skymap	:: GW skymap
		df		:: DataFrame gotten from output of Scheduler
		FOV		:: Field-of-view of the telescope
		event	:: (injections) Location of the event
	
		fr		:: Frame rate of the video (Default=4)
		'''	
		pointedTiles = self.df['Tile_Index']
		timeStamps = Time(self.df['Observation_Time']).gps.astype('int')
		previousTileIndices = np.array([], dtype='int')
		os.system('mkdir -p schedule_output')
		for times, thisTileindex, title in zip(timeStamps, pointedTiles, self.df['Observation_Time']):
			print Time(times, format='gps').iso
			print thisTileindex, previousTileIndices
			previousTileIndices = self.plotter(thisTileindex, previousTileIndices, title=title, tag=str(times))
			
		outputFileName = 'observation_schedule'
		if tag:
			outputFileName = outputFileName + '_' + tag
		command = "ffmpeg -f image2 -r " + str(fr) + " -pattern_type glob -i 'schedule_output/*.png' " + outputFileName + ".mp4"
		os.system(command)
		if not keepdata:
			os.system('rm -rf schedule_output')











