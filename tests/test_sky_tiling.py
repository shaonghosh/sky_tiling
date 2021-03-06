import numpy as np
import rankedTilesGenerator

skymap_file = 'bayestar.fits.gz'
tile_file = None

tileObj = rankedTilesGenerator.RankedTileGenerator(skymap_file, 'config.ini')
[ranked_tile_indices, ranked_tile_probs] = tileObj.getRankedTiles(resolution=256)
tileObj.plotTiles(ranked_tile_indices, ranked_tile_probs, tile_file, FOV=2.7, tileEdges=True, CI=0.99, save=True)
