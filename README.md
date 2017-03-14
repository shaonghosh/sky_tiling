# sky_tiling
The main code is the rankedTileGeneratory.py which has all the necessary methods required to find the set of tiles required to cover a given sky-loclziation. Here are the list of methods and what they do:

1.  searchedArea: Given a sky-map and the actual injection position, we can find 
                  the total searched area and the searched probability to reach 
                  to the source. The resolution of the sky-map can also be supplied.

2.  plot:         Plots the sky-map

3.  ZTF_RT:       This method generates the ranked tiles for the ZTF telescope

4.  sourceTile:   Give the ranked tiles and the actual injection position, this 
                  method finds the source tile index. 
                  
5.  detectability:  Function - given the rank of the source tile, integrtion time 
                    per tile, the total observation time, the distance to the 
                    source and limiting magnitude - time data for the telescope, 
                    this function gives the detectiability of the source




