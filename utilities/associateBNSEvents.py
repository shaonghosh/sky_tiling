import numpy as np
from astropy.io import ascii
import os




def associate(eventID, path):
        '''
        This function takes as input the event ID and the path to the ascii file and
        returns the RA, Dec and distance of the injection, and the path to the skymap.
        Make sure that the skymaps are in the same level of directory structure as the
        simulation ascii file.
        '''
        
        data = ascii.read(os.path.join(path, '2016_simulated_events.asc'))

        index = data['coinc-event-id'].data == eventID
        if np.sum(index) == 0:
        	return [np.nan, np.nan, None]

        event_RA = (data['RAdeg'].data)[index][0]
        event_Dec = (data['DEdeg'].data)[index][0]
        event_dist = (data['distance'].data)[index][0]
        event_time = (data['MJD'].data)[index][0]
        
        fitsFile = os.path.join(path, '2016_fits/') + str(eventID) + '/bayestar.fits.gz'
        return [event_time, event_RA, event_Dec, event_dist, fitsFile]