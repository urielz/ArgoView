def plot_profiles(prof,p,t,s,cyy,cmm,cdd):

    import config
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(1,facecolor='white',edgecolor='black')
    fig.set_figwidth=180
    fig.set_figheight=6

    ax = fig.add_subplot(1, 2, 1)
    ax.scatter(t,-p,edgecolor='r',color='r')
    ax.set_title(prof)
    plt.grid()
    plt.ylim([-np.ceil(np.max(p)),0 ])
    plt.xlabel('Temperature ('+config.deg+'C)')
    plt.ylabel('Pressure (dbar)')

    ax = fig.add_subplot(1, 2, 2)
    ax.scatter(s,-p,edgecolor='b',color='b')
    ax.set_title(prof)
    plt.grid()
    plt.ylim([-np.ceil(np.max(p)),0 ])
    plt.xlabel('Salinity')
    plt.ylabel('Pressure (dbar)')

    fig.set_size_inches(18.5, 10.5)
    plt.savefig(config.pdir+'%02d' %cyy+'%02d' %cmm+'%02d' %cdd+'/'+prof+'_' +'%02d' %cmm+'%02d' %cdd+'%02d' %cyy+'.png',dpi=300)
    plt.close()

    return;
