import os
from os.path import expandvars
import shutil
from mpl_toolkits.basemap import Basemap
import numpy as np
import ephem
import matplotlib.pyplot as plt
import matplotlib
import time
import logging
import tempfile
import subprocess

import utils

plt.ion()

############################################################

DECAM = 1.1 # DECam radius (deg)

############################################################

params = {
    #'backend': 'eps',
    'axes.labelsize': 16,
    #'text.fontsize': 12,           
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    #'text.usetex': True,       # ADW: Slow and no reason for tex right now
    #'font.family':'serif',
    #'font.serif':'Computer Modern Roman',
    #'figure.figsize': fig_size,
    'font.size': 12
    }
matplotlib.rcParams.update(params)

############################################################

class DECamBasemap(Basemap):

    def __init__(self, *args, **kwargs):
        super(DECamBasemap,self).__init__(self,*args,**kwargs)

    def proj(self,lon,lat):
        """ Remove points outside of projection """
        x, y = self(np.atleast_1d(lon),np.atleast_1d(lat))
        x[x > 1e29] = None
        y[y > 1e29] = None
        #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
        return x, y

    def draw_polygon(self,filename,**kwargs):
        """ Draw a polygon footprint on this Basemap instance.
        """
        defaults=dict(color='k', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        perim = np.loadtxt(filename,dtype=[('ra',float),('dec',float)])
        xy = self.proj(perim['ra'],perim['dec'])
        self.plot(*xy,**kwargs)

    def draw_des(self,**kwargs):
        """ Draw the DES footprint on this Basemap instance.
        """
        defaults=dict(color='red', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        #filename = expandvars('$MAGLITESDIR/maglites/data/round13-poly.txt')
        filename = 'data/round13-poly.txt'
        self.draw_polygon(filename,**kwargs)

    def draw_airmass(self, observatory, airmass, npts=360, **kwargs):
        defaults = dict(color='green', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)
        ra_contour = np.zeros(npts)
        dec_contour = np.zeros(npts)
        for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, npts)):
            ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
            ra_contour[ii] = np.degrees(ra_radians)
            dec_contour[ii] = np.degrees(dec_radians)
        xy = self.proj(ra_contour, dec_contour)
        self.plot(*xy, **kwargs)
         
        self.drawZenith(observatory)

    def draw_zenith(self, observatory):
        """
        Plot a to-scale representation of the focal plane size at the zenith.
        """
        defaults = dict(color='green',alpha=0.75,lw=1.5)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        # RA and Dec of zenith
        ra_zenith, dec_zenith = np.degrees(observatory.radec_of(0, '90'))
        xy = self.proj(ra_zenith, dec_zenith)
         
        self.plot(*xy,marker='+',ms=10,mew=1.5, **kwargs)
        self.tissot(ra_zenith, dec_zenith, DECAM, 100, fc='none',**kwargs)

############################################################

def drawMoon(basemap, date):
    moon = ephem.Moon()
    moon.compute(date)
    ra_moon = np.degrees(moon.ra)
    dec_moon = np.degrees(moon.dec)

    proj = safeProj(basemap, np.array([ra_moon]), np.array([dec_moon]))

    if np.isnan(proj[0]).all() or np.isnan(proj[1]).all(): return

    basemap.scatter(*proj, color='%.2f'%(0.01 * moon.phase), edgecolor='black', s=500)
    color = 'black' if moon.phase > 50. else 'white'
    plt.text(proj[0], proj[1], '%.2f'%(0.01 * moon.phase), 
             fontsize=10, ha='center', va='center', color=color)

############################################################

def safeProj(proj, lon, lat):
    """ Remove points outside of projection """
    x, y = proj(np.atleast_1d(lon),np.atleast_1d(lat))
    x[x > 1e29] = None
    y[y > 1e29] = None
    #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
    return x, y

############################################################

def drawDES(basemap, color='blue'):
    #infile = '%s/maglites/data/round13-poly.txt'%(os.environ['MAGLITESDIR'])
    infile = 'data/round13-poly.txt'
    reader_poly = open(infile)
    lines_poly = reader_poly.readlines()
    reader_poly.close()

    ra_poly = []
    dec_poly = []
    for line in lines_poly:
        if '#' in line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        ra_poly.append(float(parts[0]))
        dec_poly.append(float(parts[1]))

    proj = safeProj(basemap, ra_poly, dec_poly)
    basemap.plot(*proj, color=color, lw=3)

############################################################

def drawAirmassContour(basemap, observatory, airmass, n=360, s=50):
    #airmass = 1. / cos(90. - altitude)
    #90 - alt = arccos(1. / airmass)
    altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)

    ra_contour = np.zeros(n)
    dec_contour = np.zeros(n)
    for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, n)):
        ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
        ra_contour[ii] = np.degrees(ra_radians)
        dec_contour[ii] = np.degrees(dec_radians)
    proj = safeProj(basemap, ra_contour, dec_contour)
    basemap.plot(*proj, color='green', lw=2)

    drawZenith(basemap, observatory)
    #ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    #ra_zenith = np.degrees(ra_zenith)
    #dec_zenith = np.degrees(dec_zenith)
    #proj = safeProj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))
    #basemap.scatter(*proj, color='green', edgecolor='none', s=s)

def drawZenith(basemap, observatory):
    """
    Plot a to-scale representation of the focal plane size at the zenith.
    """
    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)
    proj = safeProj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))

    zen_kwargs = dict(color='green',alpha=0.75,lw=1.5,zorder=1000)
    basemap.plot(*proj,marker='+',ms=10,mew=1.5, **zen_kwargs)
    basemap.tissot(ra_zenith, dec_zenith, DECAM, 100, fc='none',**zen_kwargs)

############################################################

def drawMoon(basemap, date):
    moon = ephem.Moon()
    moon.compute(date)
    ra_moon = np.degrees(moon.ra)
    dec_moon = np.degrees(moon.dec)

    proj = safeProj(basemap, np.array([ra_moon]), np.array([dec_moon]))

    if np.isnan(proj[0]).all() or np.isnan(proj[1]).all(): return

    basemap.scatter(*proj, color='%.2f'%(0.01 * moon.phase), edgecolor='black', s=500)
    color = 'black' if moon.phase > 50. else 'white'
    plt.text(proj[0], proj[1], '%.2f'%(0.01 * moon.phase), 
             fontsize=10, ha='center', va='center', color=color)

############################################################

def drawTarget(basemap, ra_target, dec_target):
    proj = safeProj(basemap, np.array([ra_target]), np.array([dec_target]))

    target_kwargs = dict(color='red', marker='o', edgecolor='none', s=100, zorder=1000)
    basemap.scatter(*proj, **target_kwargs)

############################################################

def makePlot(ra, dec, date=None, name=None, figsize=(6.,6.), dpi=80, s=50, center=None, airmass=True, moon=True, des=True):
    #figsize=(10.5,8.5)
    """
    Create map in orthographic projection
    """
    if date is None: date = ephem.now()
    if type(date) != ephem.Date:
        date = ephem.Date(date)

    observatory = utils.ctio()
    observatory.date = date
    
    #fig, ax = plt.subplots(fig='ortho', figsize=FIGSIZE, dpi=DPI)
    #fig = plt.figure('ortho')
    #ax = plt.subplots(figure=fig, figsize=FIGSIZE, dpi=DPI)
    fig = plt.figure(name, figsize=figsize, dpi=dpi)

    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)

    # Zenith position
    #lon_zen = LMC_RA; lat_zen = LMC_DEC
    lon_zen = ra_zenith; lat_zen = dec_zenith

    # Create the basemap
    proj_kwargs = dict(projection='ortho', celestial=True)
    if center is None:
        lon_0, lat_0 = -lon_zen, lat_zen # Center position
    else:
        lon_0, lat_0 = center[0], center[1]

    proj_kwargs.update(lon_0=lon_0, lat_0=lat_0)
    #print proj_kwargs
    print(proj_kwargs)
    basemap = DECamBasemap(**proj_kwargs)

    parallels = np.arange(-90.,120.,30.)
    basemap.drawparallels(parallels)
    meridians = np.arange(0.,420.,60.)
    basemap.drawmeridians(meridians)

    if des: drawDES(basemap)
    if airmass: drawAirmassContour(basemap, observatory, 2., s=s)
    if moon: drawMoon(basemap, date)
    plt.title('%s UTC'%(datestring(date)))
    
    drawTarget(basemap, ra, dec)

    #return fig, ax, basemap
    return fig, basemap

############################################################

def datestring(date,precision=4): 
    """
    Convert an ephem.Date object to a string with increased precision
    
    Parameters:
    -----------
    date      : ephem.Date object
    precision : Output precision 
    
    Returns:
    --------
    datestr   : String representation of the date
    """
    """
    date = ephem.Date(date).datetime()
    datestr = date.strftime('%Y/%m/%d %H:%M:%S')
    datestr += '{:.{precision}f}'.format(date.microsecond/1.e6,
                                         precision=precision)[1:]
    """

    if precision < 0:
        msg = "Precision must be positive."
        raise Exception(msg)

    # This is a bit annoying, but works better than converting to datetime
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],precision)
    minutes = datetuple[-2]
    hours = datetuple[-3]
    minutes += int(seconds//60)
    seconds = seconds%60.
    hours += int(minutes//60)
    minutes = minutes%60

    strtuple = datetuple[:-3]+(hours,minutes,seconds)
    width = precision+2 if precision == 0 else precision+3

    #datestr = '%d/%02d/%02d %02i:%02i:%07.4f'%strtuple
    datestr = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:0{width}.{precision}f}'
    datestr = datestr.format(*strtuple,precision=precision,width=width)

    return datestr

############################################################

def nite2utc(nite, observer=None):
    import dateutil.parser
    import dateutil.tz
    import datetime

    if observer is None:
        observer = ephem.Observer()
        observer.lon = constants.LON_CTIO
        observer.lat = constants.LAT_CTIO
        observer.elevation = constants.ELEVATION_CTIO
        
    if not isinstance(nite,datetime.datetime):
        nite = dateutil.parser.parse(str(nite))

    # Move to (local) noon
    # This depends on where the user is and isn't very safe
    nite = nite.replace(hour=12,tzinfo=dateutil.tz.tzlocal())
    utc = ephem.Date(nite - nite.utcoffset())

    # Maybe something like this instead...
    #offset = int( (observer.lon/(2*np.pi)) * 24. * 60) * 60
    #tzinfo = dateutil.tz.tzoffset('OBS',offset)
    #nite = nite.replace(hour=12,tzinfo=tzinfo)
    #utc = ephem.Date(nite - nite.utcoffset())

    observer.date = utc
    observer.horizon = '-14'
    #ret = observer.next_antitransit(ephem.Sun())
    ret = observer.next_setting(ephem.Sun(), use_center=True)

    observer.horizon = '0'
    return ret

def utc2nite(utc, observer=None):
    sun = ephem.Sun()

    if observer is None:
        observer = ephem.Observer()
        observer.lon = constants.LON_CTIO
        observer.lat = constants.LAT_CTIO
        observer.elevation = constants.ELEVATION_CTIO

    observer.date = utc

    if observer.previous_setting(sun) > observer.previous_rising(sun):
        # It's night time, use the date of the previous setting
        nite = ephem.localtime(observer.previous_setting(sun))
    else:
        # It's daytime, use the next setting
        nite = ephem.localtime(observer.next_setting(sun))

    return ephem.Date(nite)

def get_nite(utc=None, observer=None):
    """Convert from a UTC date and time to the 'nite'. 

    A 'nite' is defined by the day (UTC) at noon local time in Chile
    before observing started. This follows the usual convention of
    naming the nite after the day (local time) on which it starts.

    Parameters:
    -----------
    utc : The date/time (UTC) to calculate the nite from.
    
    Returns:
    --------
    nite : An ephem.Date object containing the nite (at sunset)

    """
    if not utc: utc = ephem.now()
    return utc2nite(utc, observer)


############################################################

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

