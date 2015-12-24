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

#import py_compile
#py_compile.compile('config.py')

# choose server
server = config.server1

# check dir structure and make directories if required
if not os.path.exists(config.odir): os.system('mkdir '+config.odir)
if not os.path.exists(config.mdir): os.system('mkdir '+config.mdir)
if not os.path.exists(config.pdir): os.system('mkdir '+config.pdir)
if not os.path.exists(config.sdir): os.system('mkdir '+config.sdir)
if not os.path.exists(config.tdir): os.system('mkdir '+config.tdir)
if not os.path.exists(config.pdir2): os.system('mkdir '+config.pdir2)
if not os.path.exists(config.sdir2): os.system('mkdir '+config.sdir2)
if not os.path.exists(config.tdir2): os.system('mkdir '+config.tdir2)


# check if local Argo index file is up to date
if not os.path.exists(config.fidx): # download Argo index file
    cmd = 'ftp -V '+server+config.fidx+'.gz'
    print cmd
    print 'Downloading Argo index file. This may take some time... '
    os.system(cmd)  # get updated argo index file
    cmd = 'gunzip -f '+config.fidx+'.gz'
    os.system(cmd)
else: # check the one you have is up-to-date
    cmd = 'curl -r 0-500 '+server+'ar_index_global_prof.txt > tmp_header.txt'
    print 'Checking timestamp of Argo index file on server... '
    print cmd
    os.system(cmd)

    f = open('tmp_header.txt','r')
    for line in f:
        if 'Date of update' in line: dlineServer=line
    f.close()

    f = open(config.fidx,'r')
    for line in f:
        if 'Date of update' in line: dlineLocal=line
    f.close()

    uptodate = dlineServer == dlineLocal

    cmd = 'rm tmp_header.txt'
    os.system(cmd)

    if not uptodate:
        cmd = 'ftp -V '+server+config.fidx+'.gz'
        print 'Updating local Argo index file. This may take some time... '
        print cmd
        os.system(cmd)  # get updated argo index file
        cmd = 'gunzip -f '+config.fidx+'.gz'
        os.system(cmd)
    else:
        print 'The local Argo index file is up-to-date.'


# search for argo floats in region and time period of interest
for ii in range (config.t0,config.t1): # loop time

    cdd = date.fromordinal(ii).day
    cmm = date.fromordinal(ii).month
    cyy = date.fromordinal(ii).year

    # get new profiles metadata
    os.system('grep nc,'+str(cyy)+'%02d' %cmm+'%02d' %cdd+' '+config.fidx+' > '+config.fidx_tmp)

    # check if new profiles are inside ROI
    f = open(config.fidx_tmp,'r')
    fcnt = 0
    d = ''
    lat, lon, prid  = ([] for i in range(3))

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
        print ('No profiles available for '+str(cyy)+'/%02d' %cmm+'/%02d' %cdd+' inside the region of interest')
    else:
        print str(fcnt)+' profiles found... downloading ...'
        # download profiles on spec day in ROI
        os.chdir(config.pdir2)
        os.system('ftp -V '+d)
        os.chdir('../../')

        # make daily map with locations of new profiles
        print ('making map with locations of new profiles...')
        plot_locations(lon,lat,prid,cyy,cmm,cdd)

        if not os.path.exists(config.pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd):
            os.system('mkdir '+config.pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd)

        f = open(config.sidx, 'a')
        i = 0

        # read new profiles
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
                print ('making profile plots for float '+prof.split('_')[0][1:])
                plot_profiles(prof,p,t,s,cyy,cmm,cdd)

                if f_upd == 1: # only make trajectory and section plots if more than 1 profile for the float is available

                    # open trajectory file for this float
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','r')
                    tlat, tlon, tkm = ([] for i in range(3))
                    for line in f1:
                        tlat.append(float(line.split(' ')[1]))
                        tlon.append(float(line.split(' ')[2]))
                        tkm.append(float(line.split(' ')[3]))
                    tkm.append(sw.dist((lat[i],tlat[-1]), (lon[i],tlon[-1]), units='km')[0]+tkm[-1])
                    tlat.append(lat[i])
                    tlon.append(lon[i])
                    f1.close()
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','a')
                    f1.write('%02d' %cyy+'%02d' %cmm+'%02d' %cdd +' '+str(lat[i])+' '+str(lon[i])+' '+str(tkm[-1])+'\n')
                    f1.close()
                    i = i + 1

                    # make trajectory plot for this float
                    print ('making trajectory plot for float '+prof.split('_')[0][1:])
                    plot_trajectory(tlat,tlon,tkm,prof.split('_')[0][1:])

                    # make section plot
                    print ('making section plot for float '+prof.split('_')[0][1:])
                    fig = plt.figure(1,facecolor='white',edgecolor='black')
                    fig.set_figwidth=180
                    fig.set_figheight=6
                    ax = fig.add_subplot(2, 1, 1)
                    plot_section(prof.split('_')[0][1:])
                    plt.savefig(config.sdir+prof.split('_')[0][1:]+'.png',dpi=300)
                    plt.close()

                else: # no plots, just initialize trajectory file

                    dist_ini = 0.
                    f1 = open(config.tdir2+prof.split('_')[0][1:]+'_traj.txt','a')
                    f1.write('%02d' %cyy+'%02d' %cmm+'%02d' %cdd +' '+str(lat[i])+' '+str(lon[i])+' '+str(dist_ini)+'\n')
                    f1.close()
                    i = i + 1

            else:
                print ('No valid data in '+prof)

        f.close()

    os.system('rm '+config.fidx_tmp)
