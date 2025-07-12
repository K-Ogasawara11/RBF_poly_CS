import numpy as np
import md_cube as cube


N=80000
lon,lat,x,y,z=cube.bauer(N)

node=np.array([lon,lat,x,y,z]).T;print(node.shape)
np.savetxt('knn_bauer_N{}.txt'.format(N),node)