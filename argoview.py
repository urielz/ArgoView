# argoview.py
#
# gets and plots Argo floats profiles for a given region and time period
# initial release- dec2015
#

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import os
import datetime
import os.path
import sys
import math
from mpl_toolkits.basemap import Basemap

def geo_dist(lat1, long1, lat2, long2):

    deg2rad = math.pi/180.0 # Convert latitude and longitude to spherical coordinates in radians

    phi1 = (90.0 - lat1)*deg2rad # phi = 90 - latitude
    phi2 = (90.0 - lat2)*deg2rad

    theta1 = long1*deg2rad # theta = longitude
    theta2 = long2*deg2rad

    # Compute spherical distance from spherical coordinates
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )

    return arc

# ARGO float repositories - tries server1 first
server1  = 'ftp://usgodae.org/pub/outgoing/argo/'
server2  = 'ftp://ftp.ifremer.fr/ifremer/argo'

server = server1

# default file names
fidx     = 'ar_index_global_prof.txt'  # argo index file
fidx_tmp = 'ar_index_global_prof_tmp.txt'
sidx     = 'ar_index_sections.txt'

# directories for data (argo.*) and figures (fig.*)
mdir     = 'fig.daily_maps/'
pdir     = 'fig.daily_profiles/'
sdir     = 'fig.sections/'
tdir     = 'fig.trajectories/'
pdir2    = 'argo.profiles/'
sdir2    = 'argo.sections/'
tdir2    = 'argo.trajectories/'

# ROI boundaries latN, latS, lonW, lonE
bound = (-35,-70,-105,-35)

deg = u'\N{DEGREE SIGN}'

# check dir structure and make directories if required
if not os.path.exists(mdir): os.system('mkdir '+mdir)
if not os.path.exists(pdir): os.system('mkdir '+pdir)
if not os.path.exists(sdir): os.system('mkdir '+sdir)
if not os.path.exists(tdir): os.system('mkdir '+tdir)
if not os.path.exists(pdir2): os.system('mkdir '+pdir2)
if not os.path.exists(sdir2): os.system('mkdir '+sdir2)
if not os.path.exists(tdir2): os.system('mkdir '+tdir2)

# update argo index file
# get last file update time

if not os.path.exists(fidx): # download Argo index file
    cmd = 'ftp -V '+server+fidx+'.gz'
    print cmd
    print 'Downloading Argo index file. This may take some time... '
    os.system(cmd)  # get updated argo index file
    cmd = 'gunzip -f '+fidx+'.gz'
    os.system(cmd)
else: # check the one you have is up-to-date
    cmd = 'curl -r 0-500 '+server+'ar_index_global_prof.txt > header.txt'
    print 'Checking timestamp of Argo index file on server... '
    print cmd
    os.system(cmd)

    f = open('header.txt','r')
    for line in f:
        if 'Date of update' in line: dlineServer=line
    f.close()

    f = open(fidx,'r')
    for line in f:
        if 'Date of update' in line: dlineLocal=line
    f.close()

    uptodate = dlineServer == dlineLocal

    cmd = 'rm header.txt'
    os.system(cmd)

    if not uptodate:
        cmd = 'ftp -V '+server+fidx+'.gz'
        print 'Updating your Argo index file. This may take some time... '
        print cmd
        os.system(cmd)  # get updated argo index file
        cmd = 'gunzip -f '+fidx+'.gz'
        os.system(cmd)
    else:
        print 'Your Argo index file is up-to-date.'


# search for argo floats in ROI
for ii in range (1,31):

    d = datetime.datetime.now()
    cyy = getattr(d,'year')
    cmm = getattr(d,'month')
    cdd = getattr(d,'day')
    cmm = 11
    cdd = ii

    if os.path.exists(fidx_tmp):
        print 'current Argo file exists, execution stopped...'
#        sys.exit()

    # get new profiles metadata
    os.system('grep nc,'+str(cyy)+'%02d' %cmm+'%02d' %cdd+' '+fidx+' > '+fidx_tmp)

    f = open(fidx_tmp,'r')
    fcnt = 0
    d = ''
    lat   = []
    lon   = []
    prid  = []

    for line in f:
        clat = float(line.split(',')[2])
        clon = float(line.split(',')[3])
        cprid = str(str(line.split('/')[3]).split('.nc')[0])

        if clat < bound[0] and clat > bound[1] and clon > bound[2] and clon < bound[3]:
            lat.append(clat)
            lon.append(clon)
            prid.append(cprid)
            fcnt = fcnt + 1
            d = d+' '+server1+'dac/'+line.split(',')[0]

    f.close()

    if fcnt == 0:
        print 'no new floats profiles available'
    else:
        os.system('ftp -V '+d)
        os.system('mv *.nc ./'+pdir2)
        print str(fcnt)+' profile floats retrieved'

        # basemap plot - location of new profiles
        # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon: are the lat/lon values of the lower left and upper right corners of the map.
        m = Basemap(llcrnrlon=-140.,llcrnrlat=-70.,urcrnrlon=-40.,urcrnrlat=-20.,
                    projection='lcc',lat_1=-30.,lat_2=-60.,lon_0=-90.,
                    resolution ='l',area_thresh=1000.)

        m.bluemarble()
        parallels = np.arange(0.,-90,-5.)
        meridians = np.arange(180.,360.,15.)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

        x, y = m(lon,lat)
        m.scatter(x,y,10,marker='o',color='r')

        for label, xpt, ypt in zip(prid, x, y):
            plt.text(xpt+40000, ypt+20000, label,color='w',size='6')

        plt.title("Argo profiles - "+str(cmm)+'/'+str(cdd)+'/'+str(cyy))

        plt.savefig('./'+mdir+'/Argo_'+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd+'.png',dpi=300)
        plt.close()

        if not os.path.exists(pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd):
            os.system('mkdir '+pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd)

        f = open(sidx, 'a')

        i = 0
        for prof in prid: # read and plot individual profiles

            print prof+'.nc'

            geb = nc.Dataset('./'+pdir2+'/'+prof+'.nc')
            s = geb.variables['PSAL'][:]#.compressed()
            p = geb.variables['PRES'][:]#.compressed()
            t = geb.variables['TEMP'][:]#.compressed()

            if np.size(p)>0: # profile has valid data

                f1 = open (tdir2+prof.split('_')[0][1:]+'_traj.txt','a')
                f1.write('%02d' %cyy+'%02d' %cmm+'%02d' %cdd +' '+str(lat[i])+' '+str(lon[i])+'\n') # update trajectory file
                f1.close()
                i = i + 1

                if os.system('grep '+prof.split('_')[0][1:]+' '+sidx)!=0: # check if float is in database
                    f_upd = 0
                    f.write(prof.split('_')[0][1:]+'\n') # new float

                    with open(sdir2+prof.split('_')[0][1:]+'_t.txt','a') as f_handle:
                        np.savetxt(f_handle, p, delimiter=' ',fmt='%4d',newline='\r\n')
                        np.savetxt(f_handle, t, delimiter=' ',fmt='%3.1f',newline='\r\n')
                    with open(sdir2+prof.split('_')[0][1:]+'_s.txt','a') as f_handle:
                        np.savetxt(f_handle, p, delimiter=' ',fmt='%5d',newline='\r\n')
                        np.savetxt(f_handle, s, delimiter=' ',fmt='%3.2f',newline='\r\n')

                else: # float in database
                    f_upd = 1 # update plots

                    with open(sdir2+prof.split('_')[0][1:]+'_t.txt','a') as f_handle:
                        np.savetxt(f_handle, t, delimiter=' ',fmt='%3.1f',newline='\r\n')
                    with open(sdir2+prof.split('_')[0][1:]+'_s.txt','a') as f_handle:
                        np.savetxt(f_handle, s, delimiter=' ',fmt='%3.2f',newline='\r\n')

                fig = plt.figure(1,facecolor='white',edgecolor='black')
                fig.set_figwidth=180
                fig.set_figheight=6
                ax = fig.add_subplot(1, 2, 1)
                ax.scatter(t,-p,edgecolor='r',color='r')
                ax.set_title(prof)
                plt.grid()
                plt.ylim([-np.ceil(np.max(p)),0 ])
                plt.xlabel('Temperature ('+deg+'C)')
                plt.ylabel('Pressure (dbar)')
                ax = fig.add_subplot(1, 2, 2)
                ax.scatter(s,-p,edgecolor='b',color='b')
                ax.set_title(prof)
                plt.grid()
                plt.ylim([-np.ceil(np.max(p)),0 ])
                plt.xlabel('Salinity')
                plt.ylabel('Pressure (dbar)')
                fig.set_size_inches(18.5, 10.5)
                plt.savefig(pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd+'/'+prof+'_' +'%02d' %cmm+'%02d' %cdd+'%02d' %cyy+'.png',dpi=300)
                plt.close()

                if f_upd == 1:

                    # float trajectory

                    f1 = open (tdir2+prof.split('_')[0][1:]+'_traj.txt','r')
                    tlat = []
                    tlon = []
                    tcnt = []
                    tkm  = []
                    ii   = 0
                    for line in f1:
                        tlat.append(float(line.split(' ')[1]))
                        tlon.append(float(line.split(' ')[2]))
                        tcnt.append(ii)
                        if ii == 0:
                            tkm.append(0)
                        else:
                            if tlat[ii-1]!=tlat[ii] or tlon[ii-1]!=tlon[ii]:
                                tkm.append(6373*geo_dist(tlat[ii-1],tlon[ii-1],tlat[ii],tlon[ii])+tkm[-1])
                            else:
                                tkm.append(tkm[-1])
                        ii=ii+1
                    f1.close()

                    m = Basemap(projection='merc',llcrnrlat=max(tlat)+0.4,urcrnrlat=min(tlat)-0.4,
                    llcrnrlon=min(tlon)-0.4,urcrnrlon=max(tlon)+0.4)

                    m.bluemarble()
                    parallels = np.arange(0.,-90,-0.2)
                    meridians = np.arange(180.,360.,0.2)
                    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
                    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

                    x, y = m(tlon,tlat)
                    m.scatter(x,y,10,marker='o',color='r')
                    m.plot(x,y,'r-.')

                    lab = []
                    for j in range(1,np.size(tkm)+1):
                        lab = lab + [str(j)+' - '+str(round(tkm[j-1]))+' km']

                    for label, xpt, ypt in zip(lab, x, y):
                        plt.text(xpt+1000, ypt+500, label,color='w',size='8')

                    plt.title("Argo profile: "+str(prof.split('_')[0][1:]))

                    plt.savefig(tdir+prof.split('_')[0][1:]+'.png',dpi=300)
                    plt.close()

            else:

                print 'no valid data'

        f.close()

    os.system('rm '+fidx_tmp)
