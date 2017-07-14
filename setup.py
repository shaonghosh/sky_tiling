import os
import sys
import importlib
from utilities import createTileCenters, preComputeMap

## append this list with new packages ##
dependencies = ['argparse', 'astropy', 'healpy', 'scipy']

count = 0
for module in dependencies:
	try:
		importlib.import_module(module)
	except ImportError:
		count += 1
		print '\n--- Could not find ' + module + ' ---\n'
		continue
	print '\n*** Found package ' + module + ' ***\n'

if count > 0:
	print '\nNumber of packages not found = ' + str(count)
	raise Exception('All dependencies are NOT satisfied\n')

print '\n All dependencies are satisfied!\n'

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-P", "--path", action="store",
 				   help="Full path where the pixel-tile maps will be saved")
parser.add_argument("-W", "--work", action="store", help="Full path of the work directory")

## Only ATLAS, BlackGEM, Pan Starrs 1, ZTF precomputed files are present ##
parser.add_argument("-T", "--telescope", action="store", help="Name of telecope")
parser.add_argument("-F", "--fov", action="store", type=float, default=None, help="value of FOV, if telescope name is not among standard")
parser.add_argument("-N", "--nside", action="store", type=int, default=256, help="Target nside value")

parser.add_argument("-t", "--tilefile", action="store", help="Full path + name of tilefile if supplied by user")


## from astropy.coordinates import EarthLocation ##
## astropy.coordinates.EarthLocation.get_site_names() gives the name of valid sites ##
parser.add_argument("-S", "--site", action="store", help="Name of location of telescope")
parser.add_argument("-m", "--timemag", action="store",
				   help="Name of integration-time vs limiting-magnitude file")
parser.add_argument("-e", "--extension", action="store", help="png, pdf etc")
				   

args = parser.parse_args()

##########################################################################################

currentDir = os.getcwd()

Telescopes = ['Atlas', 'BlackGEM', 'PS1', 'ZTF'] ## List of 'standard' telescopes 


nonStandardResolution = args.nside != 256

nonStandardTelescope = args.telescope not in Telescopes

if nonStandardTelescope + nonStandardResolution:
	if args.tilefile: # check if the tile-center file is already provided by the user
		tilefile = args.tilefile
	else: # If not, then create the file, FOV is required
		if args.fov:
			tilefile = createTileCenters.createTileCenters(args.telescope, args.fov)
		else:
			print '\nFor non-standard telescopes: ' + str(Telescopes) + ' or nside =/= 256'
			print 'user must provide a valid FOV for your telescope or provide a tile-center file'
			print 'use options --fov or --tilefile when running setup.py'
			print 'Exiting...\n'
			sys.exit(1)
			
	
	### NOTE: Make sure to check for the existence of the tile-pixel map file before calling this function
	target_nside = args.nside
	if target_nside not in [64, 128, 256, 512, 1024, 2048]:
		print '\nValue of nside must be between [64, 2048] and a power of 2 '
		print 'Please provide a valid nside...'
		print 'Exiting...\n'
		sys.exit(1)
		
	precomputeFile_new = preComputeMap.preComputeMap(tilefile, args.telescope, target_nside=target_nside) # Create the tile-pixel maps 



### Preparing texts for config file ###
preComputed_64_line = 'preComputed_64 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_64.dat'
preComputed_128_line = 'preComputed_128 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_128.dat'
preComputed_256_line = 'preComputed_256 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_256.dat'
preComputed_512_line = 'preComputed_512 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_512.dat'
preComputed_1024_line = 'preComputed_1024 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_1024.dat'
preComputed_2048_line = 'preComputed_2048 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_2048.dat'


if args.tilefile:
	tileFile_line = 'tileFile = ' + args.tilefile

else:
	tileFile_line = 'tileFile = ' + currentDir + '/tile_center_files/' + args.telescope + '_tiles_indexed.dat'


filenametag_line = 'filenametag = ' + args.telescope
extension_line = 'extension = ' + args.extension

site_line = 'site = ' + args.site
time_magnitude_line = 'time_magnitude = ' + args.timemag


os.system('mkdir -p ' + args.work)
configFile = open(args.work + '/config.ini', 'w')

configFile.writelines('[pixelTileMap]\n')
configFile.writelines(preComputed_64_line + '\n')
configFile.writelines(preComputed_128_line + '\n')
configFile.writelines(preComputed_256_line + '\n')
configFile.writelines(preComputed_512_line + '\n')
configFile.writelines(preComputed_1024_line + '\n')
configFile.writelines(preComputed_2048_line + '\n\n')

configFile.writelines('[tileFiles]\n')
configFile.writelines(tileFile_line + '\n\n')

configFile.writelines('[plot]\n')
configFile.writelines(filenametag_line + '\n')
configFile.writelines(extension_line + '\n\n')

configFile.writelines('[observation]\n')
configFile.writelines(site_line + '\n')
configFile.writelines(time_magnitude_line + '\n')
configFile.writelines('trigger_time = ' + '\n\n')

configFile.close()


# if args.path is not None:
# 	os.system('cp tile_pixel_maps/*.dat ' + args.path + '/.')


binDir = currentDir + '/bin'
utils = currentDir + '/utilities'

exportText1 = 'export PYTHONPATH='+ currentDir +':${PYTHONPATH}'
exportText2 = 'export PYTHONPATH='+ binDir +':${PYTHONPATH}'
exportText3 = 'export PYTHONPATH='+ utils +':${PYTHONPATH}'
print '''\n***** sky_tiling is configured *****.
Run the following in your terminal or put it in your .bashrc'''
print exportText1
print exportText2
print exportText3
print '\n'

sourceFile = open(args.work + '/sky_tilingrc', 'w')
sourceFile.writelines(exportText1 + '\n')
sourceFile.writelines(exportText2 + '\n')
sourceFile.writelines(exportText3 + '\n')
os.system('chmod 777 ' + args.work + '/sky_tilingrc')

print 'A copy of source script is put in the work directory: ' + args.work + '/sky_tilingrc'







	


