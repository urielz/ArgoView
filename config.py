# user defined parameters for argoview.py

from datetime import date

# Initial and final dates (yyyy/mm/dd)
t0 =  date.toordinal(date(2015,11,07))
t1 =  date.toordinal(date(2015,12,21))

# Region of intereset (ROI) boundaries latN, latS, lonW, lonE
bound = (-35,-70,-105,-35)

# no of times to try to download profiles before giving up
trytimes = 2

# Argo floats global repositories
server1  = 'ftp://usgodae.org/pub/outgoing/argo/'
server2  = 'ftp://ftp.ifremer.fr/ifremer/argo'

# Argo index filename
fidx     = 'ar_index_global_prof.txt'

# Output directories for data (argo.*) and figures (fig.*)
odir = 'out/'
mdir     = odir+'fig.daily_maps/'
pdir     = odir+'fig.daily_profiles/'
sdir     = odir+'fig.sections/'
tdir     = odir+'fig.trajectories/'
pdir2    = odir+'argo.profiles/'
sdir2    = odir+'argo.sections/'
tdir2    = odir+'argo.trajectories/'

# Output filenames
fidx_tmp = odir+'ar_index_global_prof_tmp.txt'
sidx     = odir+'ar_index_sections.txt'

# GEBCO 30'' file name
gebf    = 'GEBCO_2014_2D_-105.0_-70.0_-35.0_-30.0.nc'

# Gross sanity check - max/min accepted values for variables
fillval    = 99999.
smin_check = 20.
smax_check = 40.

deg = u'\N{DEGREE SIGN}'
