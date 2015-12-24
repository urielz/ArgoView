def plot_section(profname,arg):

    import config
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

    f_p  = open(config.sdir2+profname+'_p.txt','r')
    f_xy = open(config.tdir2+profname+'_traj.txt','r')

    if arg == 0:
        f_var  = open(config.sdir2+profname+'_t.txt','r')
    if arg == 1:
        f_var  = open(config.sdir2+profname+'_s.txt','r')

    p_arr, s_arr, t_arr, lat_arr, lon_arr, dist_arr = ([] for i in range(6))

    p_cnt = 0
    for p_line in f_p: # Loop over lines and extract variables of interest
        p_line = p_line.strip().split() # removes /n (strip) and converts str to list (split)
        s_line = f_var.readline().strip().split()
        t_line = f_xy.readline().strip().split()
        lat_arr.append(float(t_line[1]))
        lon_arr.append(float(t_line[2]))
        dist_arr.append(float(t_line[3]))

        k = 0

        for i in p_line:
            if float(i)<9000:
                p_arr.append(float(i))
                s_arr.append(float(s_line[k]))
                t_arr.append(dist_arr[-1])
                k = k + 1
        if p_cnt == 0:
            np_elem = np.size(p_arr)
        p_cnt = p_cnt + 1
    f_p.close()
    f_var.close()

    p_arr = -np.asarray(p_arr)
    s_arr = np.asarray(s_arr)
    t_arr = np.asarray(t_arr)

    # define grid.
    xi = np.linspace(0,max(t_arr),p_cnt)
    yi = np.linspace(min(p_arr),0,np_elem)

    # grid the data.
    zi = griddata((t_arr, p_arr), s_arr, (xi[None,:], yi[:,None]), method='linear')
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)

    plt.colorbar() # draw colorbar
    plt.scatter(t_arr,p_arr,marker='o',c='b',s=1) # draw locations of observations
    plt.grid()
    plt.ylabel('Pressure (dbar)',fontsize=5)
    plt.xlabel('Distance (Km)',fontsize=5)

    if arg == 0:
        plt.title('Temperature section - Float: '+profname,fontsize=5)
    if arg == 1:
        plt.title('Salinity section - Float: '+profname,fontsize=5)

    return;
