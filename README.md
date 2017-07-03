# sky_tiling
The main code is the rankedTileGeneratory.py which has all the necessary methods required to find the set of tiles required to cover a given sky-loclziation. Here are the list of methods and what they do:

Create a tile object as follows:

skymap_file = 'bayestar.fits.gz' ## Name of the GW sky-localization fits file.

tileObj = rankedTilesGenerator.RankedTileGenerator(skymap_file)

To construct the ranked tiles for any given telescope, we need tile-pixel association files for that telescope. For testing purpose we include two such files in this repository, preComputed_ATLAS_pixel_indices_256.dat and preComputed_ZTF_pixel_indices_256.dat. A larger number of such files are available in NEMO cluster of UWM: /home/shaon/Tiling_Strategy/precompute_files . These are for the Atlas telescope.

"getRankedTiles" method generates the ranked tiles for any telescope whose precompted tiles-pixel maps are given.
Eample Command: [ranked_tile_indices, ranked_tile_probs] = tileObj.getRankedTiles()

These tiles could then be plotted using the method "plotTiles".
Example command: tileObj.plotTiles(ranked_tile_indices, ranked_tile_probs, FOV=None, <tile index file>, tileEdges=True)
Toggeling the tileEdges argument allows or removes the plotting of the tile boundaries. The FOV needs to be supplied in order to plot the tile boundaries.



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




