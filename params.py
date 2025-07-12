import numpy as np


# RBF setting
n=55 #stensil size
deg=5 # dgree of polydomials


# node setting
#icos
#N=2562
#N=10242 
#N=40962

#md #nomber of node
#N=1600
#N=2500
#N=6400
#N=10000
#N=25600 
#N=40000

#bauer
#N=2500
#N=6400
#N=10000
N=25600
#N=40000


#set  kind od node
#node='icos_nomod' #icosahedral grid
#node='md'  #maximum determinant
node='bauer' #spherical helix


# hyperviscosity parameter
d=2 #don't change 
k=4 # order of hyperviscosity
gamma=1e-10 # coefficient of hyperviscosity 
hep=6.6 # epsilon of hyperviscosity


# time step and output 
day=5 # number of day of integration
dt=1200 # time step in swe

day_in_sec=24*3600
total_in_sec=day*day_in_sec
hour_in_sec=60*60


rec_hour=1 # time of output in hour
steps=int(day*day_in_sec/dt)
srec0=int(hour_in_sec*rec_hour/dt)
srec=srec0
nrec=int(day*24/rec_hour)

#test settings
case=1; alpha=np.pi/4; bell='cos';lat0=0;lon0=-np.pi/2; step_in_day=120; rec_day=12;day=1200 #solid body rotation. 
# the center of bell is (lat0,lon0)=(0,-90)[deg]. time step is based on day, not on seconds.Thus, please, set step_in_day at case1. 
# rec_day determines the frequency of log creation

#case=2; alpha=np.pi/4 #Nonlinear zonal geostrophic flow (Williamosn et al. 1992)
#case=6; alpha=0 # Rossby-Haurwitz wave (Williamosn et al. 1992)

#contants
grav=9.80616
omega=7.292e-5
ea=1
rea=6.37122e6
projection='edist'
