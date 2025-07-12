# md node is ordered by kdtree.


import numpy as np
from scipy.spatial import KDTree


def load_mdorg(md):
    import numpy as np
    lad=np.loadtxt('md{}'.format(md))
    x=lad[:,0];y=lad[:,1];z=lad[:,2];wgt=lad[:,3]
    lon=np.arctan2(y,x)
    lat=np.arcsin(z)
    return lon,lat,x,y,z,wgt


N=10000
lon,lat,x,y,z,wgt=load_mdorg(N)


tree = KDTree(np.array([x,y,z]).T)
value,allindex=tree.query([x[0],y[0],z[0]],k=N)


lona=np.zeros(N)
lata=lona.copy()
xa=lona.copy()
ya=xa.copy()
za=xa.copy()
wgta=xa.copy()
for i in range(N):
    ind=int(allindex[i])
    lona[i]=lon[ind].copy()
    lata[i]=lat[ind].copy()
    xa[i]=x[ind].copy()
    ya[i]=y[ind].copy()
    za[i]=z[ind].copy()
    wgta[i]=wgt[ind].copy()

node=np.array((lona,lata,xa,ya,za,wgta));node=node.T
np.savetxt('knn_md_N{}.txt'.format(N),node)
