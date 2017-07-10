import os
import importlib

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

## from astropy.coordinates import EarthLocation ##
## astropy.coordinates.EarthLocation.get_site_names() gives the name of valid sites ##
parser.add_argument("-S", "--site", action="store", help="Name of location of telescope")
parser.add_argument("-m", "--timemag", action="store",
				   help="Name of integration-time vs limiting-magnitude file")
parser.add_argument("-e", "--extension", action="store", help="png, pdf etc")
				   

args = parser.parse_args()

currentDir = os.getcwd()

### Preparing texts for config file ###
preComputed_64 = 'preComputed_64 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_64.dat'
preComputed_128 = 'preComputed_128 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_128.dat'
preComputed_256 = 'preComputed_256 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_256.dat'
preComputed_512 = 'preComputed_512 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_512.dat'
preComputed_1024 = 'preComputed_1024 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_1024.dat'
preComputed_2048 = 'preComputed_2048 = ' + currentDir + '/tile_pixel_maps/preComputed_' + args.telescope + '_pixel_indices_2048.dat'

tileFile = 'tileFile = ' + currentDir + '/tile_center_files/' + args.telescope + '_tiles_indexed.dat'


filenametag = 'filenametag = ' + args.telescope
extension = 'extension = ' + args.extension

site = 'site = ' + args.site
time_magnitude = 'time_magnitude = ' + args.timemag


os.system('mkdir -p ' + args.work)
configFile = open(args.work + '/config.ini', 'w')

configFile.writelines('[pixelTileMap]\n')
configFile.writelines(preComputed_64 + '\n')
configFile.writelines(preComputed_128 + '\n')
configFile.writelines(preComputed_256 + '\n')
configFile.writelines(preComputed_512 + '\n')
configFile.writelines(preComputed_1024 + '\n')
configFile.writelines(preComputed_2048 + '\n\n')

configFile.writelines('[tileFiles]\n')
configFile.writelines(tileFile + '\n\n')

configFile.writelines('[plot]\n')
configFile.writelines(filenametag + '\n')
configFile.writelines(extension + '\n\n')

configFile.writelines('[observation]\n')
configFile.writelines(site + '\n')
configFile.writelines(time_magnitude + '\n')
configFile.writelines('trigger_time = ' + '\n\n')

configFile.close()





if args.path is not None:
	os.system('cp tile_pixel_maps/*.dat ' + args.path + '/.')


	
binDir = currentDir + '/bin'
exportText1 = 'export PYTHONPATH='+ currentDir +':${PYTHONPATH}'
exportText2 = 'export PYTHONPATH='+ binDir +':${PYTHONPATH}'
print '''\n***** sky_tiling is configured *****.
Run the following in your terminal or put it in your .bashrc'''
print exportText1
print exportText2
print '\n'







	


