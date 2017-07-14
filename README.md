# sky_tiling

## Requirements: ## 

0. Python 2.7
1. Argparse
2. Astropy (>1.3)
3. Healpy
4. Matplotlib
5. Numpy
6. Scipy

## Installation ##

The main ranked-tiling code is supposed to be used as a module. Thus the installation involves simply setting up the environment. This is done by the customized setup.py script. Running this code requires the following arguments:


  (*)  --work      : The name of the working directory. The configuration script will be created here.
  
  (*) --telescope : The name of the telescope. Currently only works for ATLAS, BlackGEM, Pan Starrs and ZTF. In future 
                  developments, this will work for any generic telesscopes with square FOV.
                  
  (*) --site      : The site of the observatory. This is required for the "scheduler". Note required for finding ranked-tiles.
                  However, the argument still needs to be supplied (could be just a fake site name).
                  
  (*) --extension : The file extension (png, pdf etc) of the plot of the tiles. Only required for the plotting method.
                  However, the argument still needs to be supplied (could be just a fake extension name).

Example commands:

python setup.py --work ~/RunDir/sky_tile_work --telescope BlackGEM --site None --timemag None --extension png

python setup.py --work ~/RunDir/sky_tile_work --telescope BlackGEM --site lasilla --timemag None --extension pdf

...
*NOTE: Currently the setup script has to be run from the repository itself.*
Running this should generate a config file in the work directory. The std output of the script is the environment variables that will need to be exported. Once that is done, the ranked-tiling methods will be available for use:


## Tile-Pixel maps ##

A vital ingredient for the ranked tile generation is the tile-pixel map for a telescope. This allows rapid computation of the ranked tiles. Currently the repository contains tile-pixel maps for the four aforementioned telescopes. Only resolution of nside=256 is available for all the telescopes. Higher resolution maps are larger, and we will provide them via other means (currently available in the NEMO cluster of UWM). However, the users are welcome to generate their own tile pixel maps. We describe this below:

If the user defines a telescope that is among the standard telescope name for this repository, namely Atlas, BlackGEM, PS1 and ZTF, then the ranked-tiling code will use the existing pre-computed tile-pixel maps in the utilities directory of the repository. However if a non-standard telescope name is used, then the user has two option:

Option 1.

User can specify a tile center file, which has the center of the predefined sky-grid for the tiles. The format of the tiles is as shown below (no spaces between lines). Examples of such files can be obtained from the tile_center_files directory of the repository.

ID        ra_center         dec_center

0         25.00            -89.9

1         51.43            -88.36

2         102.86           -88.36

3         154.29           -88.36

4         205.71           -88.36

The setup command that needs to be run in such a case is:

python setup.py --work ~/RunDir/sky_tile_work --telescope <some telescope> --site <some site or None> --timemag <some time magnitude file or None> --tilefile <name and full path of tile file> --extension png

This will first create the tile-pixel map file and put it in the tile_pixel_maps directory of the repository. Then create the config file in the working directory with the paths to this file in the right section.



Option 2.

If the user do not have such a file, then this file first needs to be generated first. You would need to know the field-of-view (FOV) of the telescope to generate. Note that currently the software can only handle square FOV.

To do this run the setup script with the following options:

python setup.py --work ~/RunDir/sky_tile_work --telescope <some telescope> --site <some site or None> --timemag <some time magnitude file or None> --fov <telescope FOV> --extension png

This will prompt the setup script to first generate an ad-hoc tile coordinate file that it will put in the tile_center_files directory. It will then use this file to generate the tile-pixel map files that it will put in the tile_pixel_maps directory of the repository. Finally it will create the config file and put in the work directory with the right paths.


Please note that both Option 1 and Option would require construction of the tile-pixel maps and these might take a long amount of time to run depending of the --nside value you use and the FOV of your telescope. By default (if you do not use any --nside option) the resolution used for the construction of the tile-pixel maps is nside=256. 


## Running ranked-tiling codes ##

The basic function call involves creating the ranked-tile object for a given telescope. Once the user have setup for a particular telescope, all that is required is a GW sky-localization from LIGO-Virgo. The user needs to source the environment variables that are generated by the setup script and the 'cd' into the work directory that the user defined while running setup. To test the software, the user needs to import the module rankedTilesGenerator from python and create the ranked-tile object as follows:

tileObj = rankedTilesGenerator.RankedTileGenerator('bayestar.fits.gz', 'config.ini')

where, bayestar.fits.gz is the GW sky-localization map. 

This tile object can then be used to create the ranked tiles as shown below:

[ranked_tile_indices, ranked_tile_probs] = tileObj.getRankedTiles(resolution=256)

This generates the ranked tiles and their corresponding probability values.

[Remember, to generate ranked tiles with higher resolution (unnecessary for large FOV telescopes like ZTF, Atlas or PS1), you will need to have the corresponding tile-pixel map files, which are not present in the repository. However, if you are interested to produce them, then just use a non-standard telescope name while running setup script and give the tile center file. The script will generate the tile-pixel map files].



To plot the tiles use the following method:
rankedTileData = tileObj.plotTiles(ranked_tile_indices, ranked_tile_probs, CI=0.90, FOV=47.3, tileEdges=True, save=True)


This will plot the top ranked-tiles containing the 90% localization probability. Keeping save=True will save the resulting plot in a file with the extension that the user defined while running setup script (by using the --extension option). If this save=False is used, the plot will not be saved, instead will be displayed to the user. The result of this method call is an astropy table (rankedTileData)

Rank	index	RA			Dec			Probability

1		691		24.71429	-85.93846	0.122820301363

2		733		76.14286	-85.93846	0.112929047613

3		644		127.57143	-85.93846	0.10439656724

...	

Where, the details of the 90% tiles are saved. The RA and the Dec are the location of the tile centers. The is the tile ID from the tile center file and the probability is the sum total of the probability contained within these tiles.

If tileEdges is set to False then the boundary of the tiles will not be plotted. Only the tile centers will be kept. If FOV is not given, then no tile boundary will be plotted either. 'save' allows the plot to be saved using the extension provided in the config file (generated by the setup script). CI is the confidence interval.



Other methods:


1.  searchedArea: Given a sky-map and the actual injection position, we can find 
                  the total searched area and the searched probability to reach 
                  to the source. The resolution of the sky-map can also be supplied.


2.  sourceTile:   Give the ranked tiles and the actual injection position, this 
                  method finds the source tile index. 
                  
3.  detectability:  Function - given the rank of the source tile, integrtion time 
                    per tile, the total observation time, the distance to the 
                    source and limiting magnitude - time data for the telescope, 
                    this function gives the detectiability of the source


