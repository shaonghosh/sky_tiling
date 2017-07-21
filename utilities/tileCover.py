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



def getArea(a, b, c):
    s = 0.5*(a + b + c)
    area = np.sqrt(s*(s-a)*(s-b)*(s-c))
    return area


def getTile(FOV, ra_cent, dec_cent):
    dec_down = dec_cent - 0.5*np.sqrt(FOV)
    dec_up = dec_cent + 0.5*np.sqrt(FOV)

    ra_down_left = ra_cent - 0.5*(np.sqrt(FOV)/(np.cos(dec_down*(np.pi/180.))))
    ra_down_right = ra_cent + 0.5*(np.sqrt(FOV)/(np.cos(dec_down*(np.pi/180.))))
    ra_up_left = ra_cent - 0.5*(np.sqrt(FOV)/(np.cos(dec_up*(np.pi/180.))))
    ra_up_right = ra_cent + 0.5*(np.sqrt(FOV)/(np.cos(dec_up*(np.pi/180.))))
    
    return([dec_down, dec_up, ra_down_left, ra_down_right, ra_up_left, ra_up_right])


def whichSideOfLine(x, y, x1, y1, m):
    '''
    Returns (y - y1) - m*(x - x1)
    '''
    return (y - y1) - m*(x - x1)

def slopeOfLine(x1, y1, x2, y2):
    '''
    Returns slopes corresponding to points
    '''
    return (y2 - y1)/(x2 - x1)



class GetTiles:
	def __init__(self, skymapData, FOV):
		self.FOV = FOV
		ra_map = skymapData[0]
		dec_map = skymapData[1]
		pVal = skymapData[2]
		non_zeros = pVal > 0
		self.ra = ra_map[non_zeros]
		self.dec = dec_map[non_zeros]
		self.pVal = pVal[non_zeros]
		

	def getFourTriangles(self, point_ra, point_dec):
		triangle1 = [((point_ra - self.ra_up_left)**2 + (point_dec - self.dec_up)**2)**0.5,
					((point_ra - self.ra_up_right)**2 + (point_dec - self.dec_up)**2)**0.5,
					(self.ra_up_left - self.ra_up_right)]

		triangle2 = [((point_ra - self.ra_up_left)**2 + (point_dec - self.dec_up)**2)**0.5,
					((point_ra - self.ra_down_left)**2 + (point_dec - self.dec_down)**2)**0.5,
					((self.ra_up_left - self.ra_down_left)**2 + (self.dec_up - self.dec_down)**2)**0.5]

		triangle3 = [((point_ra - self.ra_down_left)**2 + (point_dec - self.dec_down)**2)**0.5,
					((point_ra - self.ra_down_right)**2 + (point_dec - self.dec_down)**2)**0.5,
					(self.ra_down_left - self.ra_down_right)]

		triangle4 = [((point_ra - self.ra_down_right)**2 + (point_dec - self.dec_down)**2)**0.5,
					((point_ra - self.ra_up_right)**2 + (point_dec - self.dec_up)**2)**0.5,
					((self.ra_up_right - self.ra_down_right)**2 + (self.dec_up - self.dec_down)**2)**0.5]

		area1 = getArea(triangle1[0], triangle1[1], triangle1[2])
		area2 = getArea(triangle2[0], triangle2[1], triangle2[2])
		area3 = getArea(triangle3[0], triangle3[1], triangle3[2])
		area4 = getArea(triangle4[0], triangle4[1], triangle4[2])

		return area1 + area2 + area3 + area4



	def tileCover(self, tileCent, masked=None):
		'''
		masked is an array of ra and dec values whose pVals will be reset to zero from
		the sky-map. The default value of masked is None. 
		'''
		point_ra_all = self.ra
		point_dec_all = self.dec
		point_pVal_all = self.pVal
		

		[ra_cent, dec_cent] = tileCent
		[self.dec_down, self.dec_up, 
		 self.ra_down_left, self.ra_down_right, 
		 self.ra_up_left, self.ra_up_right] = getTile(self.FOV, ra_cent, dec_cent)
		 
		ra_min = np.min([self.ra_up_left, self.ra_up_right, self.ra_down_left, self.ra_down_right])
		ra_max = np.max([self.ra_up_left, self.ra_up_right, self.ra_down_left, self.ra_down_right])
		
		## Keeping part of the sky map that falls in the region in and around the tile ##
		keep = (self.dec_down <= point_dec_all)*(point_dec_all <= self.dec_up)*(ra_min <= point_ra_all)*(ra_max >= point_ra_all)
		
		point_ra_kept = point_ra_all[keep]
		point_dec_kept = point_dec_all[keep]
		point_pVal_kept = point_pVal_all[keep]

		point_ra = point_ra_kept
		point_dec = point_dec_kept
		point_pVal = point_pVal_kept

		if masked is not None:
# 			covered_earlier = np.isin(point_ra_kept, masked[:,0])*np.isin(point_dec_kept, masked[:,1])
			covered_earlier = np.in1d(point_ra_kept, masked[:,0])*np.in1d(point_dec_kept, masked[:,1])
			
			if np.sum(covered_earlier) > 0:
				point_ra = point_ra_kept[~covered_earlier]
				point_dec = point_dec_kept[~covered_earlier]
				point_pVal = point_pVal_kept[~covered_earlier]
			else:
				pass


		tile_area = self.getFourTriangles(ra_cent, dec_cent)
	
		alltriangles = self.getFourTriangles(point_ra, point_dec)
		
		if self.ra_down_left < 0:
			wrapped_point_ra = point_ra - 360.0
			alltriangles_wrapped = self.getFourTriangles(wrapped_point_ra, point_dec)
			
		if self.ra_down_right > 360:
			wrapped_point_ra = point_ra + 360.0
			alltriangles_wrapped = self.getFourTriangles(wrapped_point_ra, point_dec)


		inside = (np.round(alltriangles, 2) <= np.round(tile_area, 2))
		
		if (self.ra_down_left < 0) or (self.ra_down_right > 360):
			inside += (np.round(alltriangles_wrapped, 2) <= np.round(tile_area, 2))
				
		return [point_ra[inside], point_dec[inside], point_pVal[inside]]
		

	def tileCover_new(self, tileCent, masked=None):
		'''
		masked is an array of ra and dec values whose pVals will be reset to zero from
		the sky-map. The default value of masked is None. 
		'''
		point_ra_all = self.ra
		point_dec_all = self.dec
		point_pVal_all = self.pVal
		
		zeros = point_pVal_all == 0.0
		
		point_ra = point_ra_all[~zeros]
		point_dec = point_dec_all[~zeros]
		point_pVal = point_pVal_all[~zeros]

		if masked is not None:
			zero_these_out = np.isin(point_ra_all, masked[:,0])*np.isin(point_dec_all, masked[:,1])
			point_pVal[zero_these_out] = 0.0
		
		[ra_cent, dec_cent] = tileCent
		[self.dec_down, self.dec_up, 
		 self.ra_down_left, self.ra_down_right,
		 self.ra_up_left, self.ra_up_right] = getTile(self.FOV, ra_cent, dec_cent)
		# Only perform computation on relevant points 
		keep = (self.dec_down <= point_dec)*(point_dec <= self.dec_up)
		
		point_ra = point_ra[keep]
		point_dec = point_dec[keep]
		point_pVal = point_pVal[keep]
        
		wrapped_point_ra = point_ra - 360.0 if self.ra_down_left < 0 \
		    else point_ra + 360.0 if self.ra_down_right > 360 else point_ra
		
		slope1 = slopeOfLine(self.ra_down_left, self.dec_down,\
		     self.ra_up_left, self.dec_up)
		slope2 = slopeOfLine(self.ra_down_right, self.dec_down,\
		     self.ra_up_right, self.dec_up)
		side1 = whichSideOfLine(wrapped_point_ra, point_dec, \
		    self.ra_down_left, self.dec_down, slope1)
		side2 = whichSideOfLine(wrapped_point_ra, point_dec, \
		    self.ra_down_right, self.dec_down, slope2)
		# Pick out points which are inside the tile
		keep = ((slope1*slope2 < 0.) * (side1*side2 > 0.)) + \
		            ((slope1*slope2 > 0.) * (side1*side2 < 0.))

		return [point_ra[keep], point_dec[keep], point_pVal[keep]]








