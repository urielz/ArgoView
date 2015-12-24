def plot_trajectory(tlat,tlon,tkm,profname):

    import config
    import netCDF4 as nc
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    dlon = 0.4
    dlat = 0.4

    m = Basemap(projection='merc',llcrnrlat=max(tlat)+dlat,urcrnrlat=min(tlat)-dlat,
    llcrnrlon=min(tlon)-dlon,urcrnrlon=max(tlon)+dlon)

    geb  = nc.Dataset(config.gebf)
    topo = geb.variables['elevation'][:]
    glon = geb.variables['lon'][:]
    glat = geb.variables['lat'][:]

    # plot float trajectory
    m.bluemarble()
    parallels = np.arange(0.,-90,-0.2)
    meridians = np.arange(180.,360.,0.2)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    lon0 = np.argmin(np.abs(np.min(tlon)-dlon-glon))
    lon1 = np.argmin(np.abs(np.max(tlon)+dlon-glon))
    lat0 = np.argmin(np.abs(np.min(tlat)-dlat-glat))
    lat1 = np.argmin(np.abs(np.max(tlat)+dlat-glat))

    [LON,LAT] = np.meshgrid(glon[lon0:lon1],glat[lat0:lat1]);
    [X,Y] = m(LON,LAT)
    cmax = np.max(topo[lat0:lat1,lon0:lon1])
    cmin = np.min(topo[lat0:lat1,lon0:lon1])
    v = np.linspace(cmin,cmax,5)
    co = m.contour(X,Y,topo[lat0:lat1,lon0:lon1],v,colors='w',linestyles='solid')
    plt.clabel(co,inline=True,fmt='%1.0f',fontsize=10,colors='w')
    del X, Y

    x, y = m(tlon,tlat)
    m.scatter(x,y,10,marker='o',color='r')
    m.plot(x,y,'r-.')

    lab = []
    for j in range(1,np.size(tkm)+1):
        lab = lab + [str(j)+' - '+str(round(tkm[j-1]))+' km']

    for label, xpt, ypt in zip(lab, x, y):
        plt.text(xpt+1000, ypt+500, label,color='w',size='8')
    del x, y

    plt.title('Argo profile: '+str(profname))
    plt.savefig(config.tdir+profname+'.png',dpi=300)
    plt.close()

    return;
