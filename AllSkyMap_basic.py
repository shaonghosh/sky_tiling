from __future__ import unicode_literals

from numpy import *
import matplotlib.pyplot as pl
from matplotlib.pyplot import *
from mpl_toolkits.basemap import Basemap, pyproj

__all__ = ['AllSkyMap']

def angle_symbol(angle, round_to=1.0):
    """
    Return a string representing an angle, rounded and with a degree symbol.
    
    This is adapted from code in mpl's projections.geo module.
    """
    value = np.round(angle / round_to) * round_to
    if pl.rcParams['text.usetex'] and not pl.rcParams['text.latex.unicode']:
        return r'$%0.0f^\circ$' % value
    else:
        return '%0.0f\N{DEGREE SIGN}' % value


class AllSkyMap(Basemap):
    """
    AllSkyMap is a subclass of Basemap, specialized for handling common plotting
    tasks for celestial data.
    
    It is essentially equivalent to using Basemap with full-sphere projections
    (e.g., 'hammer' or 'moll') and the `celestial` keyword set to `True`, but
    it adds a few new methods:
    
    * label_meridians for, well, labeling meridians with their longitude values;
    
    """

    # Longitudes corresponding to east and west edges, reflecting the
    # convention that 180 deg is the eastern edge, according to basemap's 
    # underlying projections:
    east_lon = 180.
    west_lon = 180.+1.e-10

    def __init__(self, 
                       projection='hammer',
                       lat_0=0., lon_0=0.,
                       suppress_ticks=True,
                       boundinglat=None,
                       fix_aspect=True,
                       anchor=str('C'),
                       ax=None):

        if projection != 'hammer' and projection !='moll':
            raise ValueError('Only hammer and moll projections supported!')

        # Use Basemap's init, enforcing the values of many parameters that
        # aren't used or whose Basemap defaults would not be altered for all-sky
        # celestial maps.
        Basemap.__init__(self, llcrnrlon=None, llcrnrlat=None,
                       urcrnrlon=None, urcrnrlat=None,
                       llcrnrx=None, llcrnry=None,
                       urcrnrx=None, urcrnry=None,
                       width=None, height=None,
                       projection=projection, resolution=None,
                       area_thresh=None, rsphere=1.,
                       lat_ts=None,
                       lat_1=None, lat_2=None,
                       lat_0=lat_0, lon_0=lon_0,
                       suppress_ticks=suppress_ticks,
                       satellite_height=1.,
                       boundinglat=None,
                       fix_aspect=True,
                       anchor=anchor,
                       celestial=True,
                       ax=ax)

        # Keep a local ref to lon_0 for hemisphere checking.
        self._lon_0 = self.projparams['lon_0']
        self._limb = None

    def drawmapboundary(self,color='k',linewidth=1.0,fill_color=None,\
                        zorder=None,ax=None):
        """
        draw boundary around map projection region, optionally
        filling interior of region.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        line width for boundary (default 1.)
        color            color of boundary line (default black)
        fill_color       fill the map region background with this
                         color (default is no fill or fill with axis
                         background color).
        zorder           sets the zorder for filling map background
                         (default 0).
        ax               axes instance to use
                         (default None, use default axes instance).
        ==============   ====================================================

        returns matplotlib.collections.PatchCollection representing map boundary.
        """
        # Just call the base class version, but keep a copy of the limb
        # polygon for clipping.
        self._limb = Basemap.drawmapboundary(self, color=color,
            linewidth=linewidth, fill_color=fill_color, zorder=zorder, ax=ax)
        return self._limb

    def label_meridians(self, lons, fontsize=10, valign='bottom', vnudge=0,
                        halign='center', hnudge=0):
        """
        Label meridians with their longitude values in degrees.
        
        This labels meridians with negative longitude l with the value 360-l;
        for maps in celestial orientation, this means meridians to the right
        of the central meridian are labeled from 360 to 180 (left to right).
        
        `vnudge` and `hnudge` specify amounts in degress to nudge the labels
        from their default placements, vertically and horizontally.  This
        values obey the map orientation, so to nudge to the right, use a
        negative `hnudge` value.
        """
        # Run through (lon, lat) pairs, with lat=0 in each pair.
        lats = len(lons)*[0.]
        for lon,lat in zip(lons, lats):
            x, y = self(lon+hnudge, lat+vnudge)
            if lon < 0:
                lon_lbl = 360 + lon
            else:
                lon_lbl = lon
            pl.text(x, y, angle_symbol(lon_lbl), fontsize=fontsize,
                    verticalalignment=valign,
                    horizontalalignment=halign)
