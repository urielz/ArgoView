# argoview.py
#
# gets and plots Argo floats profiles for a given region and time period
#
# output:
# - Daily float locations map
# - Plot of single profiles (T and S)
# - Trajectory map for each float
# - Section map (T and S)
#
# Initial release: December 2015
#

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import datetime
from datetime import date
import os.path
import sys
import seawater as sw
import math
from mpl_toolkits.basemap import Basemap

import config
from plot_locations import plot_locations
from plot_profiles import plot_profiles
from plot_trajectory import plot_trajectory
from plot_section import plot_section
from get_index import get_index

#import py_compile
#py_compile.compile('config.py')

# choose server
server = config.server1

# Prelim 1: check dir structure and make directories if required
if not os.path.exists(config.odir): os.system('mkdir '+config.odir)
if not os.path.exists(config.mdir): os.system('mkdir '+config.mdir)
if not os.path.exists(config.pdir): os.system('mkdir '+config.pdir)
if not os.path.exists(config.sdir): os.system('mkdir '+config.sdir)
if not os.path.exists(config.tdir): os.system('mkdir '+config.tdir)
if not os.path.exists(config.pdir2): os.system('mkdir '+config.pdir2)
if not os.path.exists(config.sdir2): os.system('mkdir '+config.sdir2)
if not os.path.exists(config.tdir2): os.system('mkdir '+config.tdir2)

# Prelim 2: check timestamp of index file and decide whether to download a new one
#get_index(server)

# Main loop: search for Argo floats in region and time period of interest
for ii in range (config.t0,config.t1): # loop time

    print ('Looking for Argo profiles for '+str(date.fromordinal(ii)))

    cdd = date.fromordinal(ii).day
    cmm = date.fromordinal(ii).month
    cyy = date.fromordinal(ii).year

    # get new profiles metadata
    os.system('grep nc,'+str(cyy)+'%02d' %cmm+'%02d' %cdd+' '+config.fidx+' > '+config.fidx_tmp)

    # check if new profiles are inside ROI
    f = open(config.fidx_tmp,'r')
    fcnt = 0
    d = ''
    lat, lon, prid  = ([] for k in range(3))

    for line in f:
        clat = float(line.split(',')[2])
        clon = float(line.split(',')[3])
        cprid = str(str(line.split('/')[3]).split('.nc')[0])

        if clat < config.bound[0] and clat > config.bound[1] and clon > config.bound[2] and clon < config.bound[3]:
            lat.append(clat)
            lon.append(clon)
            prid.append(cprid)
            fcnt = fcnt + 1
            d = d+' '+server+'dac/'+line.split(',')[0]

    f.close()

    if fcnt == 0:
        print ('- No profiles available for '+str(cyy)+'/%02d' %cmm+'/%02d' %cdd+' inside the region of interest')
    else:
        print ('- Found '+str(fcnt)+' profiles... downloading ...')
        # download profiles on spec day in ROI
        os.chdir(config.pdir2)
        arfiles = d.split(' ')
        for i in range(1,np.size(arfiles)):
#            fail = os.system('ftp -V '+arfiles[i])
            if not os.path.isfile(arfiles[i].split('/')[-1]):
                fail = 1; inc = 1
                while fail and inc < config.trytimes+1:
                    fail = os.system('curl -s '+arfiles[i]+' --remote-name')
                    inc = inc + 1
                if fail:
                    print ('-- Failed to download '+arfiles[i])  #To do: decide what to do here
        os.chdir('../../')

        # make daily map with locations of new profiles
        print ('- Making map with locations of new profiles...')
        plot_locations(lon,lat,prid,cyy,cmm,cdd)

        if not os.path.exists(config.pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd):
            os.system('mkdir '+config.pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd)

        f = open(config.sidx, 'a')
        i = 0

        # read new profiles
        print ('- Load floats, make profile plots, trajectory file and sections')
        for prof in prid:

            arpr = nc.Dataset('./'+config.pdir2+'/'+prof+'.nc')
            s = arpr.variables['PSAL'][0,:]
            p = arpr.variables['PRES'][0,:]
            t = arpr.variables['TEMP'][0,:]

            p_goods = p != config.fillval

            if np.size(p_goods)>0: # profile has valid data

                # check if float is new
                if os.system('grep -q '+prof.split('_')[0][1:]+' '+config.sidx)!=0:
                    # float is not already in sections database
                    f_upd = 0 # don't make trajectory and section plots for single time!
                    f.write(prof.split('_')[0][1:]+'\n') # add float to database
                else:
                    # float in database
                    f_upd = 1 # yes, make trajectory and section plots

                # update profile files
                with open(config.sdir2+prof.split('_')[0][1:]+'_p.txt','a') as f_handle:
                    np.savetxt(f_handle, p.reshape(1,p.shape[0]), delimiter=' ',fmt='%5.1f',newline='\r\n')
                with open(config.sdir2+prof.split('_')[0][1:]+'_t.txt','a') as f_handle:
                    np.savetxt(f_handle, t.reshape(1,t.shape[0]), delimiter=' ',fmt='%3.2f',newline='\r\n')
                with open(config.sdir2+prof.split('_')[0][1:]+'_s.txt','a') as f_handle:
                    np.savetxt(f_handle, s.reshape(1,s.shape[0]), delimiter=' ',fmt='%3.2f',newline='\r\n')

                # make profile plot
#                print ('making profile plots for float '+prof.split('_')[0][1:])
                plot_profiles(prof,p,t,s,cyy,cmm,cdd)

                if f_upd == 1: # only make trajectory and section plots if more than 1 profile for the float is available

                    # open trajectory file for this float
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','r')
                    tlat, tlon, tkm = ([] for k in range(3))
                    for line in f1:
                        tlat.append(float(line.split(' ')[1]))
                        tlon.append(float(line.split(' ')[2]))
                        tkm.append(float(line.split(' ')[3]))
                    tkm.append(float(sw.dist((lat[i],tlat[-1]), (lon[i],tlon[-1]), units='km')[0])+tkm[-1])
                    tlat.append(lat[i])
                    tlon.append(lon[i])
                    f1.close()
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','a')
                    f1.write('%02d' %cyy+'%02d' %cmm+'%02d' %cdd +' '+str(lat[i])+' '+str(lon[i])+' %-5.2f' %tkm[-1]+'\n')
                    f1.close()
                    i = i + 1

                    # make trajectory plot for this float
#                    print ('making trajectory plot for float '+prof.split('_')[0][1:])
                    plot_trajectory(tlat,tlon,tkm,prof.split('_')[0][1:])

                    # make section plot
#                    print ('making section plot for float '+prof.split('_')[0][1:])
                    fig = plt.figure(1,facecolor='white',edgecolor='black')
                    fig.set_figwidth=180
                    fig.set_figheight=6
                    ax = fig.add_subplot(2, 1, 1)
                    plot_section(prof.split('_')[0][1:])
                    plt.savefig(config.sdir+prof.split('_')[0][1:]+'.png',dpi=300)
                    plt.close()

                else: # no plots, just initialize trajectory file

                    dist_ini = []
                    dist_ini.append(0.0)
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','a')
                    f1.write('%02d' %cyy+'%02d' %cmm+'%02d' %cdd +' '+str(lat[i])+' '+str(lon[i])+' '+str(dist_ini[0])+'\n')
                    f1.close()
                    i = i + 1

            else:
                print ('- No valid data in '+prof)

        f.close()

    os.system('rm '+config.fidx_tmp)
