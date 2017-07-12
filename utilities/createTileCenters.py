import numpy as np


def adjustEnd(RAs, dd, num, FOV):
    '''
    if length of RA is more than one:
    Generate a set of 'num' random numbers between zero and sqrt(FOV)/cos(decl).
    Loop over these random numbers.
    Choose a uniform set of RAs between the random value (in this iteration) and 360 degrees in intervals
    of sqrt of FOV.
    
    '''
    if len(RAs) > 1:
        randomNum = (np.sqrt(FOV)/np.cos(dd*np.pi/180.))*np.random.random(num)
        allRADiffs = np.array([])
        for ii in randomNum:
            RAs = np.arange(0+ii, 360.0, np.sqrt(FOV)/np.cos(dd*np.pi/180.))
            RAdiff = np.abs( ((360.-RAs[-1]) + RAs[0])*np.cos(dd*np.pi/180.) - (RAs[1]-RAs[0])*np.cos(dd*np.pi/180.) )
            allRADiffs = np.append(allRADiffs, RAdiff)
        II = randomNum[np.argmin(allRADiffs)]
        return np.arange(0+II, 360.0, np.sqrt(FOV)/np.cos(dd*np.pi/180.))
    else:
        return RAs


def createTileCenters(telescope, fov):
	Dec = np.arange(-90, 90, np.sqrt(fov))
	RA_vals = np.array([])
	Dec_Vals = np.array([])
	for dd in Dec:
		RAs = np.arange(0, 360.0, np.sqrt(fov)/np.cos(dd*np.pi/180.))
		newRAs = adjustEnd(RAs, dd, 10000, fov)
		Decs = dd*np.ones(len(newRAs))
		RA_vals = np.append(RA_vals, newRAs)
		Dec_Vals = np.append(Dec_Vals, Decs)
# 		gapAtEdge = ((360.0-newRAs[-1]) + newRAs[0])*np.cos(dd*np.pi/180.)
# 		avgGap = np.median(np.diff(newRAs))*np.cos(dd*np.pi/180.)
		
	allTileCents = np.vstack((np.arange(len(RA_vals)), RA_vals, Dec_Vals)).T
	tile_center_file_name = 'tile_center_files/' + telescope + '_tiles_indexed.dat'
	np.savetxt(tile_center_file_name, allTileCents, fmt='%d\t%f\t%f', header='ID\tra_center\tdec_center', comments='')
	return tile_center_file_name
		
		
		
		
		
		
		
		
		
