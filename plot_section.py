def plot_section():

#    import numpy as np
#    import matplotlib.pyplot as plt
#    from scipy.interpolate import griddata
#    import seawater as sw

    f_p = open('6901654_p.txt','r')
    f_s = open('6901654_s.txt','r')
    f_xy = open('6901654_traj.txt','r')
    p_arr=[]
    s_arr=[]
    t_arr=[]
    lat_arr=[]
    lon_arr=[]
    p_cnt = 0
    for p_line in f_p: # Loop over lines and extract variables of interest
        p_line = p_line.strip().split() # removes /n (strip) and converts str to list (split)
        s_line = f_s.readline().strip().split()
        t_line = f_xy.readline().strip().split()
        lat_arr.append(float(t_line[1]))
        lon_arr.append(float(t_line[2]))

        k = 0
        if p_cnt == 0:
            dist = 0.
        else:
            dist_prev = dist
            dist = sw.dist((lat_arr[p_cnt],lat_arr[p_cnt-1]), (lon_arr[p_cnt],lon_arr[p_cnt-1]), units='km')[0] + dist_prev
            print dist

        for i in p_line:
            if float(i)<9000:
                p_arr.append(float(i))
                s_arr.append(float(s_line[k]))
                t_arr.append(dist)
                k = k + 1
        if p_cnt == 0:
            np_elem = np.size(p_arr)
        p_cnt = p_cnt + 1
    f_p.close()
    f_s.close()

    p_arr = -np.asarray(p_arr)
    s_arr = np.asarray(s_arr)
    t_arr = np.asarray(t_arr)

    # define grid.
    xi = np.linspace(0,max(t_arr),p_cnt)
    yi = np.linspace(min(p_arr),0,np_elem)
    # grid the data.
    zi = griddata((t_arr, p_arr), s_arr, (xi[None,:], yi[:,None]), method='nearest')
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    plt.scatter(t_arr,p_arr,marker='o',c='b',s=5)

    plt.ylabel('Pressure (dbar)')
    plt.xlabel('Distance (Km)')
    plt.title("Salinity section - Float: ") #+str(prof.split('_')[0][1:]))
    #plt.savefig(sdir+prof.split('_')[0][1:]+'.png',dpi=300)
    #plt.close()

    return;
