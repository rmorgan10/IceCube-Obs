#!/usr/bin/python

"""
Docstring
"""

import sys
import numpy as np
import scipy.ndimage
#import pyfits
import healpy
import ephem
import matplotlib
import pylab

import utils
import ortho
import icecube_json

pylab.ion()

############################################################

class Event:

    def __init__(self, ra, dec, search_time=None, eventid=None, tag=None, sun_alt_threshold=-14., moon_avoidance=False, save=False):
        """
        ra and dec in degrees
        eventid is
        later in ephem Date format
        """

        self.ra = ra
        self.dec = dec
        equatorial_target = ephem.Equatorial(np.radians(self.ra), np.radians(self.dec))
        self.glon = np.degrees(ephem.Galactic(equatorial_target).lon)
        self.glat = np.degrees(ephem.Galactic(equatorial_target).lat)

        self.eventid = eventid
        self.tag = tag

        if search_time is None:
            self.search_time = ephem.now()
        else:
            self.search_time = ephem.Date(search_time)

        self.sun_alt_threshold = sun_alt_threshold

        self.observatory = utils.ctio()
        self.observatory.date = self.search_time

        if save:
            #outfile_airmass = 'pdf'
            #outfile_skymap = 'pdf'
            #outfile_ortho = 'pdf'
            #outfile_diagnostics = 'txt'
            outfile_airmass = 'png'
            outfile_skymap = 'png'
            outfile_ortho = 'png'
            outfile_diagnostics = 'txt'
        else:
            outfile_airmass = None
            outfile_skymap = None
            outfile_ortho = None
            outfile_diagnostics = None

        self.optimalTime(search_interval=1., plot=True, outfile=outfile_airmass, moon_avoidance=moon_avoidance)
        self.skymap(outfile=outfile_skymap)
        self.diagnostics(outfile=outfile_diagnostics)
        self.plotOrtho(outfile=outfile_ortho)
        if save:
            self.makeJson()

    def optimalTime(self, search_interval=1., plot=True, outfile=None, moon_avoidance=False):
        """
        Search interval in days
        """
        
        time_shift_array = np.linspace(0., search_interval, 10000)

        airmass_array = np.empty(len(time_shift_array))
        sun_alt_array = np.empty(len(time_shift_array))
        moon_alt_array = np.empty(len(time_shift_array))
        for ii, time_shift in enumerate(time_shift_array):
            date = ephem.Date(self.search_time + time_shift)
            self.observatory.date = ephem.Date(date)
        
            # Check if nighttime!
            sun_alt_array[ii] = np.degrees(ephem.Sun(self.observatory).alt)
            moon_alt_array[ii] = np.degrees(ephem.Moon(self.observatory).alt)

            ra_zenith, dec_zenith = self.observatory.radec_of(0, '90') # RA and Dec of zenith
            ra_zenith = np.degrees(ra_zenith)
            dec_zenith = np.degrees(dec_zenith)
            airmass_array[ii] = utils.getAirmass(ra_zenith, dec_zenith, self.ra, self.dec)
        
        self.observatory.date = self.search_time

        # Optimal time
        if not moon_avoidance:
            cut = (sun_alt_array < self.sun_alt_threshold)
        else:
            cut = (sun_alt_array < self.sun_alt_threshold) & (moon_alt_array < 0.)
        self.minimum_airmass = np.min(airmass_array[cut])
        optimal_time_shift = time_shift_array[cut][np.argmin(airmass_array[cut])]
        self.optimal_time = ephem.Date(self.search_time + optimal_time_shift)
            
        if plot:
            pylab.figure()
            cut = (sun_alt_array < self.sun_alt_threshold)
            
            #pylab.scatter(time_shift_array[~cut], airmass_array[~cut], c='black', s=5, edgecolor='none')
            #if np.any(cut):
            #    #pylab.scatter(time_shift_array[cut], airmass_array[cut], c=sun_alt_array[cut], edgecolor='none')
            #    #colorbar = pylab.colorbar()
            #    #colorbar.set_label('Sun Altitude (deg)')
            #    pylab.scatter(time_shift_array[cut], airmass_array[cut], c='red', edgecolor='none')
            pylab.plot(time_shift_array, airmass_array, c='black', lw=2)

            labels, n = scipy.ndimage.label(sun_alt_array > self.sun_alt_threshold)
            for index in range(1, n + 1):
                time_shift_select = time_shift_array[labels == index]
                #pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='black', alpha=0.3, zorder=-1)
                pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='red', alpha=0.15, zorder=-1)

            labels, n = scipy.ndimage.label((sun_alt_array < self.sun_alt_threshold) & (moon_alt_array > 0.))
            for index in range(1, n + 1):
                time_shift_select = time_shift_array[labels == index]
                #pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='black', alpha=0.15, zorder=-1)
                pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='yellow', alpha=0.15, zorder=-1)

            labels, n = scipy.ndimage.label((sun_alt_array < self.sun_alt_threshold) & (moon_alt_array < 0.))
            for index in range(1, n + 1):
                time_shift_select = time_shift_array[labels == index]
                #pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='black', alpha=0.15, zorder=-1)
                pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='green', alpha=0.15, zorder=-1)

            #pylab.axvline(next_rising - observatory.date, c='black', ls='--')
            #pylab.axvline(next_setting - observatory.date, c='black', ls='--')
            pylab.axvline(self.optimal_time - self.search_time, c='black', ls='--')
            pylab.axhline(self.minimum_airmass, c='black', ls='--')
            pylab.xlim(0., search_interval)
            pylab.ylim(1., 2.)
            pylab.xlabel('Date (UTC)')
            pylab.ylabel('Airmass of IceCube Field')

            xticks = [0., 0.5 * search_interval, 1. * search_interval]
            xtick_labels = [ephem.Date(self.search_time).__str__(),
                            ephem.Date(self.search_time + 0.5 * search_interval).__str__(),
                            ephem.Date(self.search_time + 1. * search_interval).__str__()]
            pylab.xticks(xticks, xtick_labels)

            if outfile:
                outfile = 'icecube_%s_airmass_%s.%s'%(self.eventid, utils.datestring(self.optimal_time), outfile)
                pylab.savefig(outfile, dpi=100)    

    def skymap(self, outfile=None):
        color_des = 'blue'

        poly_des = utils.desPoly()
        sfd = utils.openEBVMap()
        title = '%s UTC'%(self.optimal_time.__str__())
        healpy.mollview(sfd, nest=True, coord=['G','C'], min=0., max=1., xsize=3000, cmap='binary', 
                        unit='E(B-V)', title=title)
        healpy.projplot(poly_des['ra'], poly_des['dec'], lonlat=True, c=color_des, lw=2)

        healpy.projscatter(self.ra, self.dec, lonlat=True, 
                           s=85, edgecolor='white', c='red', zorder=990, clip_on=False)

        sun = ephem.Sun(self.optimal_time)
        ra_sun, dec_sun = np.degrees(sun.ra), np.degrees(sun.dec)
        healpy.projscatter(ra_sun, dec_sun, lonlat=True, 
                           s=85, edgecolor='black', c='yellow', zorder=990, clip_on=False)
    
        moon = ephem.Moon(self.optimal_time)
        ra_moon, dec_moon = np.degrees(moon.ra), np.degrees(moon.dec)
        healpy.projscatter(ra_moon, dec_moon, lonlat=True, 
                           s=85, edgecolor='black', c='%.2f'%(moon.moon_phase), zorder=990, clip_on=False)

        if outfile:
            outfile = 'output/icecube_%s_skymap_%s.%s'%(self.eventid, utils.datestring(self.optimal_time), outfile)
            pylab.savefig(outfile, dpi=100)

    def diagnostics(self, outfile=None):
        
        self.observatory.date = self.search_time

        sun = ephem.Sun(self.optimal_time)
        ra_sun, dec_sun = np.degrees(sun.ra), np.degrees(sun.dec)
        moon = ephem.Moon(self.optimal_time)
        ra_moon, dec_moon = np.degrees(moon.ra), np.degrees(moon.dec)
        sfd = utils.openEBVMap()
        nside = healpy.npix2nside(len(sfd))

        lines = ['Event',
                 '  Event ID = %s'%(self.eventid),
                 '  (ra, dec) = (%.4f, %.4f)'%(self.ra, self.dec),
                 'Date',
                 '  Now = %s (UTC)'%(ephem.now().__str__()),
                 '  Search time = %s (UTC)'%(self.search_time.__str__()),
                 '  Optimal time = %s (UTC)'%(self.optimal_time.__str__()),
                 '  Airmass at optimal time = %.2f'%(self.minimum_airmass),
                 'Sun',
                 '  Angular separation = %.2f (deg)'%(utils.angsep(self.ra, self.dec, ra_sun, dec_sun)),    
                 '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Sun()).__str__()),
                 '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Sun()).__str__()),
                 'Moon',
                 '  Illumination = %.2f'%(moon.moon_phase),
                 '  Angular separation = %.2f (deg)'%(utils.angsep(self.ra, self.dec, ra_moon, dec_moon)),
                 '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Moon()).__str__()),
                 '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Moon()).__str__()),
                 '  Next new moon = %s (UTC)'%(ephem.next_new_moon(ephem.now()).__str__()),
                 '  Next full moon = %s (UTC)'%(ephem.next_full_moon(ephem.now()).__str__()),
                 'Galactic',
                 '  (l, b) = (%.4f, %.4f)'%(self.glon, self.glat),
                 '  E(B-V) = %.2f'%(sfd[utils.angToPix(nside, self.glon, self.glat)])]
        for line in lines:
            #print line
            print(line)
        
        if outfile:
            outfile = 'output/icecube_%s_diagnostics_%s.%s'%(self.eventid, utils.datestring(self.optimal_time), outfile)
            writer = open(outfile, 'w')
            for line in lines:
                writer.write(line + '\n')
            writer.close()
        
        ###
        """
        print 'Event'
        print '  Event ID = %s'%(self.eventid)
        print '  (ra, dec) = (%.2f, %.2f)'%(self.ra, self.dec)

        print 'Date'
        print '  Now = %s (UTC)'%(ephem.now().__str__())
        print '  Search time = %s (UTC)'%(self.search_time.__str__())
        print '  Optimal time = %s (UTC)'%(self.optimal_time.__str__())
        print '  Airmass at optimal time = %.2f'%(self.minimum_airmass)

        print 'Sun'
        sun = ephem.Sun(self.optimal_time)
        ra_sun, dec_sun = np.degrees(sun.ra), np.degrees(sun.dec)
        print '  Angular separation = %.2f (deg)'%(utils.angsep(self.ra, self.dec, ra_sun, dec_sun))    
        print '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Sun()).__str__())
        print '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Sun()).__str__())

        print 'Moon'
        moon = ephem.Moon(self.optimal_time)
        ra_moon, dec_moon = np.degrees(moon.ra), np.degrees(moon.dec)
        print '  Illumination = %.2f'%(moon.moon_phase)
        print '  Angular separation = %.2f (deg)'%(utils.angsep(self.ra, self.dec, ra_moon, dec_moon))
        print '  Next rising = %s (UTC)'%(self.observatory.next_rising(ephem.Moon()).__str__())
        print '  Next setting = %s (UTC)'%(self.observatory.next_setting(ephem.Moon()).__str__())
        print '  Next new moon = %s (UTC)'%(ephem.next_new_moon(ephem.now()).__str__())
        print '  Next full moon = %s (UTC)'%(ephem.next_full_moon(ephem.now()).__str__())
    
        print 'Galactic'
        sfd = utils.openEBVMap()
        nside = healpy.npix2nside(len(sfd))
        print '  (l, b) = (%.2f, %.2f)'%(self.glon, self.glat)
        print '  E(B-V) = %.2f'%(sfd[utils.angToPix(nside, self.glon, self.glat)])
        """

    def plotOrtho(self, outfile=None):
        ortho.makePlot(self.ra, self.dec, date=self.optimal_time)
        if outfile:
            outfile = 'icecube_%s_ortho_%s.%s'%(self.eventid, utils.datestring(self.optimal_time), outfile)
            pylab.savefig(outfile)

    def makeJson(self):
        outfile = 'output/icecube_%s_%s.json'%(self.eventid,
                                        utils.datestring(self.optimal_time))
        #self.optimal_time.__str__().replace('/','-').replace(' ','-').replace(':','-'))

        sispi = icecube_json.makeJson(self.ra, self.dec, self.eventid, self.optimal_time)
        icecube_json.writeJson(outfile, sispi)

############################################################

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ra', dest='ra_target', type=float, default=None,
                        help='RA (deg)')
    parser.add_argument('--dec', dest='dec_target', type=float, default=None,
                        help='Dec (deg)')
    parser.add_argument('--eventid', type=int,
                        help='Event ID number')
    parser.add_argument('--time', default=None,
                        help='Search time in format YYYY/MM/DD HH:MM:SS')
    parser.add_argument('--tag',  default=None,
                        help='Label used for plotting and outfiles')
    parser.add_argument('--moon', action='store_true', default=False, 
                        help='Moon avoidance')
    parser.add_argument('--alt', dest='sun_alt_threshold', type=float, default=-14.,
                        help='Sun altitude threshold (deg)')
    args = parser.parse_args()
    
    #print 'Args:\n', args, '\n'
    print('Args:\n', args, '\n')
    if args.ra_target is None or args.dec_target is None or args.eventid is None:
        parser.print_help()
        sys.exit()

    event = Event(args.ra_target, args.dec_target, eventid=args.eventid,
                  search_time=args.time, tag=args.tag,
                  sun_alt_threshold=args.sun_alt_threshold, moon_avoidance=args.moon, save=True)

    raw_input('WAIT')
    sys.exit()

    """
    diagnostics(ra_target, dec_target)
    #skymap(ra_target, dec_target, outfile=tag + '_skymap.pdf')
    #timeline(ra_target, dec_target, outfile=tag + '_airmass.pdf')
    skymap(ra_target, dec_target, outfile=None)
    timeline(ra_target, dec_target, outfile=None)
    """

############################################################
"""
def timeline(ra_target, dec_target, sun_alt_threshold=-14., outfile=None):
    
    observatory = ctio()
    time_shift_array = np.linspace(0., 1., 10000)

    airmass_array = []
    sun_alt_array = []
    for time_shift in time_shift_array:
        date = ephem.Date(ephem.now() + time_shift)
        observatory.date = ephem.Date(date)
        
        # Check if nighttime!
        sun_alt_array.append(np.degrees(ephem.Sun(observatory).alt))

        ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
        ra_zenith = np.degrees(ra_zenith)
        dec_zenith = np.degrees(dec_zenith)
        airmass_array.append(getAirmass(ra_zenith, dec_zenith, ra_target, dec_target))

    airmass_array = np.array(airmass_array)
    sun_alt_array = np.array(sun_alt_array)

    # Optimal time
    cut = (sun_alt_array < sun_alt_threshold)
    optimal_airmass = np.min(airmass_array[cut])
    optimal_time_shift = time_shift_array[cut][np.argmin(airmass_array[cut])]
    optimal_time = ephem.Date(ephem.now() + optimal_time_shift)
    print 'Airmass'
    print '  Minimum airmass = %.2f'%(optimal_airmass)
    print '  Optimal time = %s'%(optimal_time.__str__())

    # Next rise and set
    observatory.date = ephem.now()
    observatory.horizon = np.radians(sun_alt_threshold)
    next_rising = observatory.next_rising(ephem.Sun())
    next_setting = observatory.next_setting(ephem.Sun())

    pylab.figure()
    cut = (sun_alt_array < sun_alt_threshold)
    pylab.scatter(time_shift_array[~cut], airmass_array[~cut], c='black', s=5, edgecolor='none')
    if np.any(cut):
        #pylab.scatter(time_shift_array[cut], airmass_array[cut], c=sun_alt_array[cut], edgecolor='none')
        #colorbar = pylab.colorbar()
        #colorbar.set_label('Sun Altitude (deg)')
        pylab.scatter(time_shift_array[cut], airmass_array[cut], c='red', edgecolor='none')

    labels, n = scipy.ndimage.label(sun_alt_array > sun_alt_threshold)
    for index in range(1, n + 1):
        time_shift_select = time_shift_array[labels == index]
        pylab.axvspan(np.min(time_shift_select), np.max(time_shift_select), color='black', alpha=0.15, zorder=-1)

    #pylab.axvline(next_rising - observatory.date, c='black', ls='--')
    #pylab.axvline(next_setting - observatory.date, c='black', ls='--')
    pylab.axvline(optimal_time - observatory.date, c='black', ls='--')
    pylab.axhline(optimal_airmass, c='black', ls='--')
    pylab.xlim(0., 1.)
    pylab.ylim(1., 2.)
    pylab.xlabel('Date (UTC)')
    pylab.ylabel('Airmass of IceCube Field')

    xticks = [0., 0.5, 1.]
    xtick_labels = [ephem.Date(ephem.now() + 0.0).__str__(),
                    ephem.Date(ephem.now() + 0.5).__str__(),
                    ephem.Date(ephem.now() + 1.0).__str__()]
    print xtick_labels
    pylab.xticks(xticks, xtick_labels)

    if outfile:
        pylab.savefig(outfile, dpi=100)
"""
############################################################
