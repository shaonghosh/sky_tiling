##################################
# Author: Shaon Ghosh            #
# Radboud University             #
# email: shaon@astro.ru.nl       #
##################################

import numpy as np
import healpy as hp

class cover:
    def __init__(self, pVal, ra, dec, tileBound, CI):
        self.p = pVal
        self.ra = ra
        self.dec = dec
        self.tileBound = tileBound
        self.CI = CI


    def contourCovering(self):
        RA = self.ra
        Dec = self.dec
        pVal = self.p
        pValSum = np.cumsum(pVal)
        include = pValSum < self.CI*np.sum(pVal)
        if np.sum(include) < len(pValSum): include[np.sum(include)] = True
        RA = RA[include]
        Dec = Dec[include]
        pValinclude = pVal[include]
              
        coverTiles = []
        pValCover = []
        for tile in self.tileBound:
            intile = np.all([RA > tile[0], RA < tile[1], Dec > tile[2], Dec < tile[3]], axis=0)
            if np.sum(intile) > 0:
                coverTiles.append(tile)
                pValCover.append(np.sum(pValinclude[intile]))
        
        tileorder = np.argsort(-np.array(pValCover))
        coverTiles = np.array(coverTiles)[tileorder]
        return coverTiles
        
    def rankedTiling(self):
        RA = self.ra
        Dec = self.dec
        pVal = self.p
        yy = []
        for pixel in self.tileBound:
            yy.append( np.all([RA > pixel[0], RA < pixel[1], Dec > pixel[2], Dec < pixel[3]], axis=0) )
            
        yy = np.array(yy)
        Sum = np.sum(yy, axis=0)
        tVal_normalized = []
        pixelinTile = []
        ValTiles = []
        
        for ii in range(0, len(pVal)):
            if Sum[ii] != 0: pVal_normalized.append(pVal[ii]/Sum[ii])
            if Sum[ii] == 0: pVal_normalized.append(pVal[ii])  ## Sanity check for incomplete grid
        tVal_normalized = np.array(tVal_normalized)
        tMatrix = np.array([np.multiply(tVal_normalized, yy[ii,:]) for ii in range(0, len(pixelBound))])

        check = np.sum(tMatrix, axis=1) > 0
        RankedTiles = self.tileBound[check]
        for intile in yy[check]:
            ValTiles.append(np.sum(tVal_normalized[intile]))
            pixelinTile.append(tVal_normalized[intile])

        pixelinTile = np.array(pixelinTile)
        tileorder = np.argsort(-np.array(pValTiles))

        ValTiles = np.array(ValTiles)[tileorder]
        pixelinTile = pixelinTile[tileorder]
        tValSum = np.cumsum(ValTiles)
        include = tValSum < self.CI*np.sum(ValTiles)
        del tValSum
        include[np.sum(include)] = True 
        pixelinTile_EB = pixelinTile[include]









