def plot_locations(lon,lat,prid,cyy,cmm,cdd):

    import config
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

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

    plt.savefig('./'+config.mdir+'Argo_'+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd+'.png',dpi=300)
    plt.close()

    return;
