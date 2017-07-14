import placeTile
import healpy as hp
import numpy as np
import tileCover
import cPickle as pickle


skymapFile = '/Users/ghosh4/Downloads/2016_fits/417956/bayestar.fits.gz'
skymapOrig = hp.read_map(skymapFile, verbose=False)
skymap = hp.ud_grade(skymapOrig, 128, power=-2)
npix = len(skymap)
nside = hp.npix2nside(npix)
theta, phi = hp.pix2ang(nside, np.arange(0, npix))
ra_map = np.rad2deg(phi)
dec_map = np.rad2deg(0.5*np.pi - theta)
pVal = skymap[np.arange(0, npix)]
order = np.argsort(-pVal)
ra_map = ra_map[order]
dec_map = dec_map[order]
pVal = pVal[order]
skymapData = [ra_map, dec_map, pVal]
FOV=47.929036686

include = np.cumsum(pVal) < 0.9
include[np.sum(include)] = True

ra_90 = ra_map[include]
dec_90 = dec_map[include]
pVal_90 = pVal[include]


File = open('samples3Tile.pickle', 'rb')
output = pickle.load(File)
ra_samples_1 = output[:,0]
dec_samples_1 = output[:,1]
ra_samples_2 = output[:,2]
dec_samples_2 = output[:,3]
ra_samples_3 = output[:,4]
dec_samples_3 = output[:,5]

numtiles=2

obj = placeTile.PlaceTile(skymapData, FOV, numtiles=numtiles)
# result = obj.localizeTC(ra_samples_1=ra_samples_1, dec_samples_1=dec_samples_1, ra_samples_2=ra_samples_2, dec_samples_2=dec_samples_2, ra_samples_3=ra_samples_3, dec_samples_3=dec_samples_3)
# x = obj.getSamples()
# print x
# print np.shape(x)
# for ii in range(0, numtiles):
# 	print x[:,ii*2]
# 	print np.max(x[:,ii*2]), np.min(x[:,ii*2])
# 	print x[:,(ii*2)+1]
# 	print np.max(x[:,(ii*2)+1]), np.min(x[:,(ii*2)+1])
# 	print '*****************'

result = obj.localizeTC()

print result[0]
print '****************'
print result[1]

plotObj = placeTile.PlotterClass(result[0], result[1], FOV=FOV)
plotObj.plotter(ra_90, dec_90, 'TwoTileCase_withOverlapCondition_' + str(numtiles) + 'Tiles.pdf')

masked_points=np.array([]).reshape(0, 2)
tileObj = tileCover.GetTiles(skymapData, FOV)
for ii in range(0, len(result[0])):
	tileCent = [result[0][ii], result[1][ii]]

	thisPointing = tileObj.tileCover(tileCent, masked_points)
	points_to_be_masked = np.vstack((thisPointing[0], thisPointing[1])).T
	masked_points = np.vstack((masked_points, points_to_be_masked))
	print 'Probability of this tile = ' + str(np.sum(thisPointing[2]))



# secondPointing = tileObj.tileCover(tileCent_2, points)
# firstTwoPointings = [np.append(firstPointing[0], secondPointing[0]), np.append(firstPointing[1], secondPointing[1])]
# points = np.vstack((firstTwoPointings[0], firstTwoPointings[1])).T
# thirdPointing = tileObj.tileCover(tileCent_3, points)
# 
# 
# 
# 
# print '****************'
# print np.sum(firstPointing[2])
# print np.sum(secondPointing[2])
# print np.sum(thirdPointing[2])








