import numpy as np
import matplotlib.pyplot as plt
import md_cube as cube
import matplotlib.tri as tri
from matplotlib import cm, colors
from params import N,n,node,case

notc=2
time=24*5

h=np.loadtxt('case{}_h_N{}.txt'.format(notc,N))
cu=np.loadtxt('case{}_cu_N{}.txt'.format(notc,N))
cv=np.loadtxt('case{}_cv_N{}.txt'.format(notc,N))

if node=='md':
    lon,lat,x,y,z,wgt=cube.knn_md(N)

elif node=='icos_nomod':
    lon,lat,x,y,z=cube.knn_icos_nomod(N)
    wgt=np.loadtxt('../nodes/knn_icos_nomodify_weight{}.txt'.format(N))
   
triang = tri.Triangulation(lon, lat)
fig=plt.figure(figsize=(4,4))
ax=fig.add_subplot(111)

if case==2:
    ht=np.loadtxt('case{}_init_h_N{}.txt'.format(notc,N))
    tcf=ax.tricontourf(triang,h-ht,levels=11,cmap='coolwarm')
    plt.colorbar(tcf,ax=ax)
    plt.savefig("case{}_h_N{}n{}_{}hour.png".format(notc,N,n,time))
    plt.show()
    
    h0=wgt@abs(ht)
    h1=wgt@abs(ht-h)
    l1=h1/h0
    h0h0=(wgt@(ht**2))**(1/2)
    dh=(wgt@((h-ht)**2))**(1/2)
    l2=dh/h0h0
    linf=np.amax(abs(ht-h))/np.amax(abs(ht))

    print('l1={:e}'.format(l1))
    print('l2={:e}'.format(l2))
    print('linf={:e}'.format(linf))
