import numpy as np
from numpy import sin,cos,tan,pi
import matplotlib.pyplot as plt
import scipy
import md_cube as cube
import md_dmatrix as dmat
import md_ordering as order
import initial_condition as init
from params import N,n,dt,steps,k,gamma,node,projection,deg,hep,omega,case,ea,rea,alpha,d
from params import srec,nrec,srec0
import os
from tqdm import tqdm
import md_integration as integration
import matplotlib.tri as tri
import time


#constant
delta=-1/(N**k)*gamma
a=ea/np.sqrt(3)


print('N{},n{},dt{},steps{}'.format(N,n,dt,steps))
print('{} node'.format(node))
print('{} projection'.format(projection))

# loading nodes
if node=='md':
    lon,lat,x,y,z,wgt=cube.knn_md(N)

elif node=='icos_nomod':
    lon,lat,x,y,z=cube.knn_icos_nomod(N)
    wgt=np.loadtxt('../nodes/knn_icos_nomodify_weight{}.txt'.format(N))

elif node =='bauer':
    lon,lat,x,y,z=cube.bauer(N)
    wgt=np.loadtxt('../nodes/bauer_weight_N{}.txt'.format(N))
   

# node partitioning

idim1,idim2,idim3,idim4,idim5,idim6=cube.partition_xyz(x,y,z,lat,lon,N)
x1=x[idim1];y1=y[idim1];z1=z[idim1]
x2=x[idim2];y2=y[idim2];z2=z[idim2]
x3=x[idim3];y3=y[idim3];z3=z[idim3]
x4=x[idim4];y4=y[idim4];z4=z[idim4]
x5=x[idim5];y5=y[idim5];z5=z[idim5]
x6=x[idim6];y6=y[idim6];z6=z[idim6]


#get index and order index
cindex1=order.convert_stensil_index(x,y,z,x1,y1,z1,idim1,1,n)
cindex2=order.convert_stensil_index(x,y,z,x2,y2,z2,idim2,2,n)
cindex3=order.convert_stensil_index(x,y,z,x3,y3,z3,idim3,3,n)
cindex4=order.convert_stensil_index(x,y,z,x4,y4,z4,idim4,4,n)
cindex5=order.convert_stensil_index(x,y,z,x5,y5,z5,idim5,5,n)
cindex6=order.convert_stensil_index(x,y,z,x6,y6,z6,idim6,6,n)

cidim1=order.ordering_idim(x1,y1,z1,idim1,1)
cidim2=order.ordering_idim(x2,y2,z2,idim2,2)
cidim3=order.ordering_idim(x3,y3,z3,idim3,3)
cidim4=order.ordering_idim(x4,y4,z4,idim4,4)
cidim5=order.ordering_idim(x5,y5,z5,idim5,5)
cidim6=order.ordering_idim(x6,y6,z6,idim6,6)

x1=x[cindex1];y1=y[cindex1];z1=z[cindex1]
x2=x[cindex2];y2=y[cindex2];z2=z[cindex2]
x3=x[cindex3];y3=y[cindex3];z3=z[cindex3]
x4=x[cindex4];y4=y[cindex4];z4=z[cindex4]
x5=x[cindex5];y5=y[cindex5];z5=z[cindex5]
x6=x[cindex6];y6=y[cindex6];z6=z[cindex6]


lon1=lon[cidim1];lat1=lat[cidim1]
lon2=lon[cidim2];lat2=lat[cidim2]
lon3=lon[cidim3];lat3=lat[cidim3]
lon4=lon[cidim4];lat4=lat[cidim4]
lon5=lon[cidim5];lat5=lat[cidim5]
lon6=lon[cidim6];lat6=lat[cidim6]

slon1=lon[cindex1];slat1=lat[cindex1]
slon2=lon[cindex2];slat2=lat[cindex2]
slon3=lon[cindex3];slat3=lat[cindex3]
slon4=lon[cindex4];slat4=lat[cindex4]
slon5=lon[cindex5];slat5=lat[cindex5]
slon6=lon[cindex6];slat6=lat[cindex6]


#find interaciton index
overlaparea=order.calc_overlaparea(cindex1,cindex2,cindex3,cindex4,cindex5,cindex6)
over2inL1,over2inR1,c2innerL1,c2bR1,b2overlap1,in2bR1,in2cR1=order.calc_initeraction_index(overlaparea,cidim1,cindex1)
over2inL2,over2inR2,c2innerL2,c2bR2,b2overlap2,in2bR2,in2cR2=order.calc_initeraction_index(overlaparea,cidim2,cindex2)
over2inL3,over2inR3,c2innerL3,c2bR3,b2overlap3,in2bR3,in2cR3=order.calc_initeraction_index(overlaparea,cidim3,cindex3)
over2inL4,over2inR4,c2innerL4,c2bR4,b2overlap4,in2bR4,in2cR4=order.calc_initeraction_index(overlaparea,cidim4,cindex4)
over2inL5,over2inR5,c2innerL5,c2bR5,b2overlap5,in2bR5,in2cR5=order.calc_initeraction_index(overlaparea,cidim5,cindex5)
over2inL6,over2inR6,c2innerL6,c2bR6,b2overlap6,in2bR6,in2cR6=order.calc_initeraction_index(overlaparea,cidim6,cindex6)

uni1,uni2,uni3,uni4,uni5,uni6=order.calc_overlaparea2(cindex1,cidim1,cindex2,cidim2,cindex3,cidim3,cindex4,cidim4,cindex5,cidim5,cindex6,cidim6)
over2inL1,over2inR1=order.calc_initeraction_indexv2(overlaparea,cidim1,cindex1,uni1)
over2inL2,over2inR2=order.calc_initeraction_indexv2(overlaparea,cidim2,cindex2,uni2)
over2inL3,over2inR3=order.calc_initeraction_indexv2(overlaparea,cidim3,cindex3,uni3)
over2inL4,over2inR4=order.calc_initeraction_indexv2(overlaparea,cidim4,cindex4,uni4)
over2inL5,over2inR5=order.calc_initeraction_indexv2(overlaparea,cidim5,cindex5,uni5)
over2inL6,over2inR6=order.calc_initeraction_indexv2(overlaparea,cidim6,cindex6,uni6)

# preparations of transformation for overlaparea
olon=lon[overlaparea]
olat=lat[overlaparea]

o2inveclat1=olat[over2inR1];o2inveclat2=olat[over2inR2];o2inveclat3=olat[over2inR3]
o2inveclat4=olat[over2inR4];o2inveclat5=olat[over2inR5];o2inveclat6=olat[over2inR6]

o2inveclon1=olon[over2inR1];o2inveclon2=olon[over2inR2];o2inveclon3=olon[over2inR3]
o2inveclon4=olon[over2inR4];o2inveclon5=olon[over2inR5];o2inveclon6=olon[over2inR6]

c2overlat1=lat1[c2bR1];c2overlat2=lat2[c2bR2];c2overlat3=lat3[c2bR3]
c2overlat4=lat4[c2bR4];c2overlat5=lat5[c2bR5];c2overlat6=lat6[c2bR6]

c2overlon1=lon1[c2bR1];c2overlon2=lon2[c2bR2];c2overlon3=lon3[c2bR3]
c2overlon4=lon4[c2bR4];c2overlon5=lon5[c2bR5];c2overlon6=lon6[c2bR6]

#set differencial matrix

if os.path.isfile('../Dmatrix/Dx1_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx1_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy1_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL1_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
    index=np.load('../Dmatrix/index1_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x1,y1,z1,cidim1,1,n,deg,projection)
    
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim1.size,n)
    np.save('../Dmatrix/Dx1_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy1_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/HFL1_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/index1_{}N{}n{}deg{}'.format(node,N,n,deg),index)
vsize=x1.size
SDx1,SDy1=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim1,N,n,vsize,c2innerL1)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim1.size,c2innerL1)
hiv1=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x1.size,x1.size))


#P2
if os.path.isfile('../Dmatrix/Dx2_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx2_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy2_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index=np.load('../Dmatrix/index2_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL2_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x2,y2,z2,cidim2,2,n,deg,projection)
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim2.size,n)
    np.save('../Dmatrix/Dx2_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy2_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index2_{}N{}n{}deg{}'.format(node,N,n,deg),index)
    np.save('../Dmatrix/HFL2_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
vsize=x2.size
SDx2,SDy2=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim2,N,n,vsize,c2innerL2)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim2.size,c2innerL2)
hiv2=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x2.size,x2.size))


#P3
if os.path.isfile('../Dmatrix/Dx3_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx3_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy3_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index=np.load('../Dmatrix/index3_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL3_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x3,y3,z3,cidim3,3,n,deg,projection)
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim3.size,n)
    np.save('../Dmatrix/HFL3_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/Dx3_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy3_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index3_{}N{}n{}deg{}'.format(node,N,n,deg),index)
vsize=x3.size
SDx3,SDy3=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim3,N,n,vsize,c2innerL3)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim3.size,c2innerL3)
hiv3=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x3.size,x3.size))



#P4
if os.path.isfile('../Dmatrix/Dx4_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx4_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy4_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index=np.load('../Dmatrix/index4_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL4_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x4,y4,z4,cidim4,4,n,deg,projection)
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim4.size,n)
    np.save('../Dmatrix/HFL4_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/Dx4_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy4_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index4_{}N{}n{}deg{}'.format(node,N,n,deg),index)
vsize=x4.size
SDx4,SDy4=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim4,N,n,vsize,c2innerL4)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim4.size,c2innerL4)
hiv4=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x4.size,x4.size))


#P5
if os.path.isfile('../Dmatrix/Dx5_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx5_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy5_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index=np.load('../Dmatrix/index5_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL5_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x5,y5,z5,cidim5,5,n,deg,projection)
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim5.size,n)
    np.save('../Dmatrix/HFL5_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/Dx5_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy5_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index5_{}N{}n{}deg{}'.format(node,N,n,deg),index)
vsize=x5.size
SDx5,SDy5=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim5,N,n,vsize,c2innerL5)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim5.size,c2innerL5)
hiv5=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x5.size,x5.size))


#P6
if os.path.isfile('../Dmatrix/Dx6_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx6_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy6_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index=np.load('../Dmatrix/index6_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    H=np.load('../Dmatrix/HFL6_{}N{}n{}ep{:.1f}K{}.npy'.format(node,N,n,hep,k))
else:
    Dx,Dy,r,index,R=dmat.calc_Dmatrix_phs(x,y,z,x6,y6,z6,cidim6,6,n,deg,projection)
    H=dmat.calc_hiv_FH2011(d,k,hep,R,cidim6.size,n)
    np.save('../Dmatrix/HFL6_{}N{}n{}ep{:.1f}K{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/Dx6_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy6_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index6_{}N{}n{}deg{}'.format(node,N,n,deg),index)
vsize=x6.size
SDx6,SDy6=dmat.convert_sparse_c2inner(index,Dx/rea,Dy/rea,cidim6,N,n,vsize,c2innerL6)
row,col,vals=dmat.make_row_col_c2input(index,H*delta,n,cidim6.size,c2innerL6)
hiv6=scipy.sparse.csr_matrix((vals,(row,col)),shape = (x6.size,x6.size))


col1=2*omega*(-cos(slon1)*cos(slat1)*sin(alpha)+sin(slat1)*cos(alpha))
col2=2*omega*(-cos(slon2)*cos(slat2)*sin(alpha)+sin(slat2)*cos(alpha))
col3=2*omega*(-cos(slon3)*cos(slat3)*sin(alpha)+sin(slat3)*cos(alpha))
col4=2*omega*(-cos(slon4)*cos(slat4)*sin(alpha)+sin(slat4)*cos(alpha))
col5=2*omega*(-cos(slon5)*cos(slat5)*sin(alpha)+sin(slat5)*cos(alpha))
col6=2*omega*(-cos(slon6)*cos(slat6)*sin(alpha)+sin(slat6)*cos(alpha))


# jakobian of metric tensor
gi1=cos(slat1)**3*cos(slon1)**3/a**2

slon2b=slon2-pi/2
gi2=cos(slat2)**3*cos(slon2b)**3/a**2


slon3b=slon3.copy()
for i in range(slon3.size):
    if slon3[i]<0:
        slon3b[i]=slon3b[i]+pi*2
slon3b-=pi
gi3=cos(slat3)**3*cos(slon3b)**3/a**2

slon4b=slon4+pi/2
gi4=cos(slat4)**3*cos(slon4b)**3/a**2

gi5=sin(slat5)**3/a**2
gi6=np.sqrt(sin(slat6)**6/a**4)

#sin,cos,tan of lat and lon
sinlon1=np.sin(slon1);coslon1=cos(slon1);tanlon1=tan(slon1)
sinlon2=np.sin(slon2b);coslon2=cos(slon2b);tanlon2=tan(slon2b)
sinlon3=np.sin(slon3b);coslon3=cos(slon3b);tanlon3=tan(slon3b)
sinlon4=np.sin(slon4b);coslon4=cos(slon4b);tanlon4=tan(slon4b)
sinlon5=sin(slon5);coslon5=cos(slon5)
sinlon6=sin(slon6);coslon6=cos(slon6)

sinlat1=sin(slat1);coslat1=cos(slat1);tanlat1=tan(slat1)
sinlat2=sin(slat2);coslat2=cos(slat2);tanlat2=tan(slat2)
sinlat3=sin(slat3);coslat3=cos(slat3);tanlat3=tan(slat3)
sinlat4=sin(slat4);coslat4=cos(slat4);tanlat4=tan(slat4)
sinlat5=sin(slat5);coslat5=cos(slat5)
sinlat6=sin(slat6);coslat6=cos(slat6)

sin_c2overlat1=sin(c2overlat1);cos_c2overlat1=cos(c2overlat1)
sin_c2overlat2=sin(c2overlat2);cos_c2overlat2=cos(c2overlat2)
sin_c2overlat3=sin(c2overlat3);cos_c2overlat3=cos(c2overlat3)
sin_c2overlat4=sin(c2overlat4);cos_c2overlat4=cos(c2overlat4)
sin_c2overlat5=sin(c2overlat5);cos_c2overlat5=cos(c2overlat5)
sin_c2overlat6=sin(c2overlat6);cos_c2overlat6=cos(c2overlat6)

sin_c2overlon5=sin(c2overlon5);cos_c2overlon5=cos(c2overlon5)
sin_c2overlon6=sin(c2overlon6);cos_c2overlon6=cos(c2overlon6)


sin_c2overlon1=sin(c2overlon1);cos_c2overlon1=cos(c2overlon1)

c2overlonb=c2overlon2.copy()-pi/2
sin_c2overlon2=sin(c2overlonb);cos_c2overlon2=cos(c2overlonb)

c2overlonb=c2overlon3.copy()
for i in range(c2overlon3.size):
    if c2overlon3[i]<0:
        c2overlonb[i]+=pi*2
c2overlonb-=pi
sin_c2overlon3=sin(c2overlonb);cos_c2overlon3=cos(c2overlonb)

c2overlonb=c2overlon4.copy()+pi/2
sin_c2overlon4=sin(c2overlonb);cos_c2overlon4=cos(c2overlonb)

#over2difvec
sin_o2inveclat1=sin(o2inveclat1);cos_o2inveclat1=cos(o2inveclat1);tan_o2inveclat1=np.tan(o2inveclat1)
sin_o2inveclat2=sin(o2inveclat2);cos_o2inveclat2=cos(o2inveclat2);tan_o2inveclat2=np.tan(o2inveclat2)
sin_o2inveclat3=sin(o2inveclat3);cos_o2inveclat3=cos(o2inveclat3);tan_o2inveclat3=np.tan(o2inveclat3)
sin_o2inveclat4=sin(o2inveclat4);cos_o2inveclat4=cos(o2inveclat4);tan_o2inveclat4=np.tan(o2inveclat4)
sin_o2inveclat5=sin(o2inveclat5);cos_o2inveclat5=cos(o2inveclat5)
sin_o2inveclat6=sin(o2inveclat6);cos_o2inveclat6=cos(o2inveclat6)

sin_o2inveclon5=sin(o2inveclon5);cos_o2inveclon5=cos(o2inveclon5)
sin_o2inveclon6=sin(o2inveclon6);cos_o2inveclon6=cos(o2inveclon6)


sin_o2inveclon1=sin(o2inveclon1);cos_o2inveclon1=cos(o2inveclon1);tan_o2inveclon1=tan(o2inveclon1)

o2inveclonb=o2inveclon2.copy()-pi/2
sin_o2inveclon2=sin(o2inveclonb);cos_o2inveclon2=cos(o2inveclonb);tan_o2inveclon2=tan(o2inveclonb)

o2inveclonb=o2inveclon3.copy()
for i in range(o2inveclon3.size):
    if o2inveclon3[i]<0:
        o2inveclonb[i]+=pi*2
o2inveclonb-=pi
sin_o2inveclon3=sin(o2inveclonb);cos_o2inveclon3=cos(o2inveclonb);tan_o2inveclon3=tan(o2inveclonb)

o2inveclonb=o2inveclon4.copy()+pi/2
sin_o2inveclon4=sin(o2inveclonb);cos_o2inveclon4=cos(o2inveclonb);tan_o2inveclon4=tan(o2inveclonb)


# get initial conditions
if case==2:
    u,v,h=init.case2(lon,lat)

elif case==6:
    u,v,h=init.case6(lon,lat)

log_h=np.zeros((nrec+1,N))
log_cu=np.zeros((nrec+1,N))
log_cv=np.zeros((nrec+1,N))
log_h[0]=h.copy()

#input initial conditions
hi1=np.zeros(cindex1.size);hi2=np.zeros(cindex2.size);hi3=np.zeros(cindex3.size)
hi4=np.zeros(cindex4.size);hi5=np.zeros(cindex5.size);hi6=np.zeros(cindex6.size)

hi1[c2innerL1]=h[cidim1];hi2[c2innerL2]=h[cidim2];hi3[c2innerL3]=h[cidim3]
hi4[c2innerL4]=h[cidim4];hi5[c2innerL5]=h[cidim5];hi6[c2innerL6]=h[cidim6]

ho=h[np.int64(overlaparea)]

ho0=ho.copy()
us=u[overlaparea]
vs=v[overlaparea]

# calc contravariant

u1,v1=cube.ca2cs_eqv2(u[cindex1],v[cindex1],coslat1,coslon1,tanlat1,tanlon1)
u2,v2=cube.ca2cs_eqv2(u[cindex2],v[cindex2],coslat2,coslon2,tanlat2,tanlon2)
u3,v3=cube.ca2cs_eqv2(u[cindex3],v[cindex3],coslat3,coslon3,tanlat3,tanlon3)
u4,v4=cube.ca2cs_eqv2(u[cindex4],v[cindex4],coslat4,coslon4,tanlat4,tanlon4)
u5,v5=cube.ca2cs_npv2(u[cindex5],v[cindex5],coslat5,coslon5,sinlat5,sinlon5)
u6,v6=cube.ca2cs_spv2(u[cindex6],v[cindex6],coslat6,coslon6,sinlat6,sinlon6)


#calc covariant
cui1,cvi1=cube.calc_covariant_eqv2(u1,v1,coslat1,coslon1,sinlat1,sinlon1)
cui2,cvi2=cube.calc_covariant_eqv2(u2,v2,coslat2,coslon2,sinlat2,sinlon2)
cui3,cvi3=cube.calc_covariant_eqv2(u3,v3,coslat3,coslon3,sinlat3,sinlon3)
cui4,cvi4=cube.calc_covariant_eqv2(u4,v4,coslat4,coslon4,sinlat4,sinlon4)
cui5,cvi5=cube.calc_covariant_npv2(u5,v5,coslat5,coslon5,sinlat5,sinlon5)
cui6,cvi6=cube.calc_covariant_spv2(u6,v6,coslat6,coslon6,sinlat6,sinlon6)

acu=np.zeros(N)
acv=acu.copy()
acu[cidim1]=cui1[c2innerL1];acv[cidim1]=cvi1[c2innerL1]
acu[cidim2]=cui2[c2innerL2];acv[cidim2]=cvi2[c2innerL2]
acu[cidim3]=cui3[c2innerL3];acv[cidim3]=cvi3[c2innerL3]
acu[cidim4]=cui4[c2innerL4];acv[cidim4]=cvi4[c2innerL4]
acu[cidim5]=cui5[c2innerL5];acv[cidim5]=cvi5[c2innerL5]
acu[cidim6]=cui6[c2innerL6];acv[cidim6]=cvi6[c2innerL6]

log_cu[0,:]=acu
log_cv[0,:]=acv



#RK4
du1=np.zeros(cidim1.size);dv1=np.zeros(cidim1.size);dh1=np.zeros(cidim1.size)
du2=np.zeros(cidim2.size);dv2=np.zeros(cidim2.size);dh2=np.zeros(cidim2.size)
du3=np.zeros(cidim3.size);dv3=np.zeros(cidim3.size);dh3=np.zeros(cidim3.size)
du4=np.zeros(cidim4.size);dv4=np.zeros(cidim4.size);dh4=np.zeros(cidim4.size)
du5=np.zeros(cidim5.size);dv5=np.zeros(cidim5.size);dh5=np.zeros(cidim5.size)
du6=np.zeros(cidim6.size);dv6=np.zeros(cidim6.size);dh6=np.zeros(cidim6.size)

count=1  
start=time.time()
for i in tqdm(range(steps)):
#for i in range(steps):

    h1c=hi1.copy();h2c=hi2.copy();h3c=hi3.copy();h4c=hi4.copy();h5c=hi5.copy();h6c=hi6.copy()
    cu1c=cui1.copy();cu2c=cui2.copy();cu3c=cui3.copy();cu4c=cui4.copy();cu5c=cui5.copy();cu6c=cui6.copy()
    cv1c=cvi1.copy();cv2c=cvi2.copy();cv3c=cvi3.copy();cv4c=cvi4.copy();cv5c=cvi5.copy();cv6c=cvi6.copy()
    
    # calc contravariant components
    ui1,vi1=cube.co2contra_eqv2(cu1c,cv1c,coslat1,coslon1,sinlat1,sinlon1)
    ui2,vi2=cube.co2contra_eqv2(cu2c,cv2c,coslat2,coslon2,sinlat2,sinlon2)
    ui3,vi3=cube.co2contra_eqv2(cu3c,cv3c,coslat3,coslon3,sinlat3,sinlon3)
    ui4,vi4=cube.co2contra_eqv2(cu4c,cv4c,coslat4,coslon4,sinlat4,sinlon4)
    ui5,vi5=cube.co2contra_npv2(cu5c,cv5c,coslat5,coslon5,sinlat5,sinlon5)
    ui6,vi6=cube.co2contra_spv2(cu6c,cv6c,coslat6,coslon6,sinlat6,sinlon6)
   
    #overlapvector to difvec
    hi1[over2inL1]=ho[over2inR1];hi2[over2inL2]=ho[over2inR2];hi3[over2inL3]=ho[over2inR3]
    hi4[over2inL4]=ho[over2inR4];hi5[over2inL5]=ho[over2inR5];hi6[over2inL6]=ho[over2inR6]
    
     
    us1=us[over2inR1];us2=us[over2inR2];us3=us[over2inR3];us4=us[over2inR4];us5=us[over2inR5];us6=us[over2inR6]
    vs1=vs[over2inR1];vs2=vs[over2inR2];vs3=vs[over2inR3];vs4=vs[over2inR4];vs5=vs[over2inR5];vs6=vs[over2inR6]
   
    u1b,v1b=cube.ca2cs_eqv2(us1,vs1,cos_o2inveclat1,cos_o2inveclon1,tan_o2inveclat1,tan_o2inveclon1);ui1[over2inL1]=u1b;vi1[over2inL1]=v1b
    u2b,v2b=cube.ca2cs_eqv2(us2,vs2,cos_o2inveclat2,cos_o2inveclon2,tan_o2inveclat2,tan_o2inveclon2);ui2[over2inL2]=u2b;vi2[over2inL2]=v2b
    u3b,v3b=cube.ca2cs_eqv2(us3,vs3,cos_o2inveclat3,cos_o2inveclon3,tan_o2inveclat3,tan_o2inveclon3);ui3[over2inL3]=u3b;vi3[over2inL3]=v3b
    u4b,v4b=cube.ca2cs_eqv2(us4,vs4,cos_o2inveclat4,cos_o2inveclon4,tan_o2inveclat4,tan_o2inveclon4);ui4[over2inL4]=u4b;vi4[over2inL4]=v4b
    u5b,v5b=cube.ca2cs_npv2(us5,vs5,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);ui5[over2inL5]=u5b;vi5[over2inL5]=v5b
    u6b,v6b=cube.ca2cs_spv2(us6,vs6,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);ui6[over2inL6]=u6b;vi6[over2inL6]=v6b
    
    cu1b,cv1b=cube.calc_covariant_eqv2(u1b,v1b,cos_o2inveclat1,cos_o2inveclon1,sin_o2inveclat1,sin_o2inveclon1);cui1[over2inL1]=cu1b;cvi1[over2inL1]=cv1b
    cu2b,cv2b=cube.calc_covariant_eqv2(u2b,v2b,cos_o2inveclat2,cos_o2inveclon2,sin_o2inveclat2,sin_o2inveclon2);cui2[over2inL2]=cu2b;cvi2[over2inL2]=cv2b
    cu3b,cv3b=cube.calc_covariant_eqv2(u3b,v3b,cos_o2inveclat3,cos_o2inveclon3,sin_o2inveclat3,sin_o2inveclon3);cui3[over2inL3]=cu3b;cvi3[over2inL3]=cv3b
    cu4b,cv4b=cube.calc_covariant_eqv2(u4b,v4b,cos_o2inveclat4,cos_o2inveclon4,sin_o2inveclat4,sin_o2inveclon4);cui4[over2inL4]=cu4b;cvi4[over2inL4]=cv4b
    cu5b,cv5b=cube.calc_covariant_npv2(u5b,v5b,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);cui5[over2inL5]=cu5b;cvi5[over2inL5]=cv5b
    cu6b,cv6b=cube.calc_covariant_spv2(u6b,v6b,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);cui6[over2inL6]=cu6b;cvi6[over2inL6]=cv6b    

    
    du1,dv1,dh1=integration.calc_tendencyv2(SDx1,SDy1,hiv1,cui1,cvi1,ui1,vi1,hi1,gi1,col1,dt)
    du2,dv2,dh2=integration.calc_tendencyv2(SDx2,SDy2,hiv2,cui2,cvi2,ui2,vi2,hi2,gi2,col2,dt)
    du3,dv3,dh3=integration.calc_tendencyv2(SDx3,SDy3,hiv3,cui3,cvi3,ui3,vi3,hi3,gi3,col3,dt)
    du4,dv4,dh4=integration.calc_tendencyv2(SDx4,SDy4,hiv4,cui4,cvi4,ui4,vi4,hi4,gi4,col4,dt)
    du5,dv5,dh5=integration.calc_tendencyv2(SDx5,SDy5,hiv5,cui5,cvi5,ui5,vi5,hi5,gi5,col5,dt)
    du6,dv6,dh6=integration.calc_tendencyv2(SDx6,SDy6,hiv6,cui6,cvi6,ui6,vi6,hi6,gi6,col6,dt)


    hi1+=dh1/6;hi2+=dh2/6;hi3+=dh3/6
    hi4+=dh4/6;hi5+=dh5/6;hi6+=dh6/6

    cui1+=du1/6;cui2+=du2/6;cui3+=du3/6
    cui4+=du4/6;cui5+=du5/6;cui6+=du6/6

    cvi1+=dv1/6;cvi2+=dv2/6;cvi3+=dv3/6
    cvi4+=dv4/6;cvi5+=dv5/6;cvi6+=dv6/6
    

    
#K2
    
    
    cut1=cu1c+du1/2;cut2=cu2c+du2/2;cut3=cu3c+du3/2
    cut4=cu4c+du4/2;cut5=cu5c+du5/2;cut6=cu6c+du6/2

    cvt1=cv1c+dv1/2;cvt2=cv2c+dv2/2;cvt3=cv3c+dv3/2
    cvt4=cv4c+dv4/2;cvt5=cv5c+dv5/2;cvt6=cv6c+dv6/2

    ht1=h1c+dh1/2;ht2=h2c+dh2/2;ht3=h3c+dh3/2
    ht4=h4c+dh4/2;ht5=h5c+dh5/2;ht6=h6c+dh6/2
    
    #input vals from planevals to overlapvector
    hb1=ht1[in2bR1];ho[b2overlap1]=hb1.copy();hb2=ht2[in2bR2];ho[b2overlap2]=hb2.copy()
    hb3=ht3[in2bR3];ho[b2overlap3]=hb3.copy();hb4=ht4[in2bR4];ho[b2overlap4]=hb4.copy()
    hb5=ht5[in2bR5];ho[b2overlap5]=hb5.copy();hb6=ht6[in2bR6];ho[b2overlap6]=hb6.copy()

    #overlap to difvec

    # calc contravariant components
    ui1,vi1=cube.co2contra_eqv2(cut1,cvt1,coslat1,coslon1,sinlat1,sinlon1)
    ui2,vi2=cube.co2contra_eqv2(cut2,cvt2,coslat2,coslon2,sinlat2,sinlon2)
    ui3,vi3=cube.co2contra_eqv2(cut3,cvt3,coslat3,coslon3,sinlat3,sinlon3)
    ui4,vi4=cube.co2contra_eqv2(cut4,cvt4,coslat4,coslon4,sinlat4,sinlon4)
    ui5,vi5=cube.co2contra_npv2(cut5,cvt5,coslat5,coslon5,sinlat5,sinlon5)
    ui6,vi6=cube.co2contra_spv2(cut6,cvt6,coslat6,coslon6,sinlat6,sinlon6)
    u1b=ui1[in2bR1];u2b=ui2[in2bR2];u3b=ui3[in2bR3];u4b=ui4[in2bR4];u5b=ui5[in2bR5];u6b=ui6[in2bR6]
    v1b=vi1[in2bR1];v2b=vi2[in2bR2];v3b=vi3[in2bR3];v4b=vi4[in2bR4];v5b=vi5[in2bR5];v6b=vi6[in2bR6]
      
    
    #transform  vector form contravariant to sphere
    us1,vs1=cube.cs2ca_eqv2(u1b,v1b,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
    us2,vs2=cube.cs2ca_eqv2(u2b,v2b,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
    us3,vs3=cube.cs2ca_eqv2(u3b,v3b,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
    us4,vs4=cube.cs2ca_eqv2(u4b,v4b,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
    us5,vs5=cube.cs2ca_npv2(u5b,v5b,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
    us6,vs6=cube.cs2ca_spv2(u6b,v6b,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)
    
    #input windvector in ovelapvector
    us[b2overlap1]=us1.copy();us[b2overlap2]=us2.copy();us[b2overlap3]=us3.copy();us[b2overlap4]=us4.copy();us[b2overlap5]=us5.copy();us[b2overlap6]=us6.copy()
    vs[b2overlap1]=vs1.copy();vs[b2overlap2]=vs2.copy();vs[b2overlap3]=vs3.copy();vs[b2overlap4]=vs4.copy();vs[b2overlap5]=vs5.copy();vs[b2overlap6]=vs6.copy()
    
    
    # transform vals in overlapvector
    ht1[over2inL1]=ho[over2inR1];ht2[over2inL2]=ho[over2inR2];ht3[over2inL3]=ho[over2inR3]
    ht4[over2inL4]=ho[over2inR4];ht5[over2inL5]=ho[over2inR5];ht6[over2inL6]=ho[over2inR6]

    # input  vals in overlaparea to difvec
    us1=us[over2inR1];us2=us[over2inR2];us3=us[over2inR3];us4=us[over2inR4];us5=us[over2inR5];us6=us[over2inR6]
    vs1=vs[over2inR1];vs2=vs[over2inR2];vs3=vs[over2inR3];vs4=vs[over2inR4];vs5=vs[over2inR5];vs6=vs[over2inR6]
    u1b,v1b=cube.ca2cs_eqv2(us1,vs1,cos_o2inveclat1,cos_o2inveclon1,tan_o2inveclat1,tan_o2inveclon1);ui1[over2inL1]=u1b;vi1[over2inL1]=v1b
    u2b,v2b=cube.ca2cs_eqv2(us2,vs2,cos_o2inveclat2,cos_o2inveclon2,tan_o2inveclat2,tan_o2inveclon2);ui2[over2inL2]=u2b;vi2[over2inL2]=v2b
    u3b,v3b=cube.ca2cs_eqv2(us3,vs3,cos_o2inveclat3,cos_o2inveclon3,tan_o2inveclat3,tan_o2inveclon3);ui3[over2inL3]=u3b;vi3[over2inL3]=v3b
    u4b,v4b=cube.ca2cs_eqv2(us4,vs4,cos_o2inveclat4,cos_o2inveclon4,tan_o2inveclat4,tan_o2inveclon4);ui4[over2inL4]=u4b;vi4[over2inL4]=v4b
    u5b,v5b=cube.ca2cs_npv2(us5,vs5,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);ui5[over2inL5]=u5b;vi5[over2inL5]=v5b
    u6b,v6b=cube.ca2cs_spv2(us6,vs6,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);ui6[over2inL6]=u6b;vi6[over2inL6]=v6b

    cu1b,cv1b=cube.calc_covariant_eqv2(u1b,v1b,cos_o2inveclat1,cos_o2inveclon1,sin_o2inveclat1,sin_o2inveclon1);cut1[over2inL1]=cu1b;cvt1[over2inL1]=cv1b
    cu2b,cv2b=cube.calc_covariant_eqv2(u2b,v2b,cos_o2inveclat2,cos_o2inveclon2,sin_o2inveclat2,sin_o2inveclon2);cut2[over2inL2]=cu2b;cvt2[over2inL2]=cv2b
    cu3b,cv3b=cube.calc_covariant_eqv2(u3b,v3b,cos_o2inveclat3,cos_o2inveclon3,sin_o2inveclat3,sin_o2inveclon3);cut3[over2inL3]=cu3b;cvt3[over2inL3]=cv3b
    cu4b,cv4b=cube.calc_covariant_eqv2(u4b,v4b,cos_o2inveclat4,cos_o2inveclon4,sin_o2inveclat4,sin_o2inveclon4);cut4[over2inL4]=cu4b;cvt4[over2inL4]=cv4b
    cu5b,cv5b=cube.calc_covariant_npv2(u5b,v5b,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);cut5[over2inL5]=cu5b;cvt5[over2inL5]=cv5b
    cu6b,cv6b=cube.calc_covariant_spv2(u6b,v6b,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);cut6[over2inL6]=cu6b;cvt6[over2inL6]=cv6b
    
    # calc_tendency
    du1,dv1,dh1=integration.calc_tendencyv2(SDx1,SDy1,hiv1,cut1,cvt1,ui1,vi1,ht1,gi1,col1,dt)
    du2,dv2,dh2=integration.calc_tendencyv2(SDx2,SDy2,hiv2,cut2,cvt2,ui2,vi2,ht2,gi2,col2,dt)
    du3,dv3,dh3=integration.calc_tendencyv2(SDx3,SDy3,hiv3,cut3,cvt3,ui3,vi3,ht3,gi3,col3,dt)
    du4,dv4,dh4=integration.calc_tendencyv2(SDx4,SDy4,hiv4,cut4,cvt4,ui4,vi4,ht4,gi4,col4,dt)
    du5,dv5,dh5=integration.calc_tendencyv2(SDx5,SDy5,hiv5,cut5,cvt5,ui5,vi5,ht5,gi5,col5,dt)
    du6,dv6,dh6=integration.calc_tendencyv2(SDx6,SDy6,hiv6,cut6,cvt6,ui6,vi6,ht6,gi6,col6,dt)
   
    
    
    hi1+=dh1/3;hi3+=dh3/3
    hi4+=dh4/3;hi2+=dh2/3
    hi5+=dh5/3;hi6+=dh6/3
    cui1+=du1/3;cui2+=du2/3;cui3+=du3/3
    cui4+=du4/3;cui5+=du5/3;cui6+=du6/3
    cvi1+=dv1/3;cvi2+=dv2/3;cvi3+=dv3/3
    cvi4+=dv4/3;cvi5+=dv5/3;cvi6+=dv6/3
    
    
    
    #K3
    cut1=cu1c+du1/2;cut2=cu2c+du2/2;cut3=cu3c+du3/2
    cut4=cu4c+du4/2;cut5=cu5c+du5/2;cut6=cu6c+du6/2

    cvt1=cv1c+dv1/2;cvt2=cv2c+dv2/2;cvt3=cv3c+dv3/2
    cvt4=cv4c+dv4/2;cvt5=cv5c+dv5/2;cvt6=cv6c+dv6/2

    ht1=h1c+dh1/2;ht2=h2c+dh2/2;ht3=h3c+dh3/2
    ht4=h4c+dh4/2;ht5=h5c+dh5/2;ht6=h6c+dh6/2
    
    #input vals from planevals to overlapvector
    hb1=ht1[in2bR1];ho[b2overlap1]=hb1.copy();hb2=ht2[in2bR2];ho[b2overlap2]=hb2.copy()
    hb3=ht3[in2bR3];ho[b2overlap3]=hb3.copy();hb4=ht4[in2bR4];ho[b2overlap4]=hb4.copy()
    hb5=ht5[in2bR5];ho[b2overlap5]=hb5.copy();hb6=ht6[in2bR6];ho[b2overlap6]=hb6.copy()

    
     
    # calc contravariant components
    ui1,vi1=cube.co2contra_eqv2(cut1,cvt1,coslat1,coslon1,sinlat1,sinlon1)
    ui2,vi2=cube.co2contra_eqv2(cut2,cvt2,coslat2,coslon2,sinlat2,sinlon2)
    ui3,vi3=cube.co2contra_eqv2(cut3,cvt3,coslat3,coslon3,sinlat3,sinlon3)
    ui4,vi4=cube.co2contra_eqv2(cut4,cvt4,coslat4,coslon4,sinlat4,sinlon4)
    ui5,vi5=cube.co2contra_npv2(cut5,cvt5,coslat5,coslon5,sinlat5,sinlon5)
    ui6,vi6=cube.co2contra_spv2(cut6,cvt6,coslat6,coslon6,sinlat6,sinlon6)
    u1b=ui1[in2bR1];u2b=ui2[in2bR2];u3b=ui3[in2bR3];u4b=ui4[in2bR4];u5b=ui5[in2bR5];u6b=ui6[in2bR6]
    v1b=vi1[in2bR1];v2b=vi2[in2bR2];v3b=vi3[in2bR3];v4b=vi4[in2bR4];v5b=vi5[in2bR5];v6b=vi6[in2bR6]
      
    
    #transform  vector form contravariant to sphere
    us1,vs1=cube.cs2ca_eqv2(u1b,v1b,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
    us2,vs2=cube.cs2ca_eqv2(u2b,v2b,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
    us3,vs3=cube.cs2ca_eqv2(u3b,v3b,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
    us4,vs4=cube.cs2ca_eqv2(u4b,v4b,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
    us5,vs5=cube.cs2ca_npv2(u5b,v5b,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
    us6,vs6=cube.cs2ca_spv2(u6b,v6b,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)
    

    us[b2overlap1]=us1.copy();us[b2overlap2]=us2.copy();us[b2overlap3]=us3.copy();us[b2overlap4]=us4.copy();us[b2overlap5]=us5.copy();us[b2overlap6]=us6.copy()
    vs[b2overlap1]=vs1.copy();vs[b2overlap2]=vs2.copy();vs[b2overlap3]=vs3.copy();vs[b2overlap4]=vs4.copy();vs[b2overlap5]=vs5.copy();vs[b2overlap6]=vs6.copy()
    
    
    # overlaparea to input vector
    ht1[over2inL1]=ho[over2inR1];ht2[over2inL2]=ho[over2inR2];ht3[over2inL3]=ho[over2inR3]
    ht4[over2inL4]=ho[over2inR4];ht5[over2inL5]=ho[over2inR5];ht6[over2inL6]=ho[over2inR6]

    us1=us[over2inR1];us2=us[over2inR2];us3=us[over2inR3];us4=us[over2inR4];us5=us[over2inR5];us6=us[over2inR6]
    vs1=vs[over2inR1];vs2=vs[over2inR2];vs3=vs[over2inR3];vs4=vs[over2inR4];vs5=vs[over2inR5];vs6=vs[over2inR6]
    u1b,v1b=cube.ca2cs_eqv2(us1,vs1,cos_o2inveclat1,cos_o2inveclon1,tan_o2inveclat1,tan_o2inveclon1);ui1[over2inL1]=u1b;vi1[over2inL1]=v1b
    u2b,v2b=cube.ca2cs_eqv2(us2,vs2,cos_o2inveclat2,cos_o2inveclon2,tan_o2inveclat2,tan_o2inveclon2);ui2[over2inL2]=u2b;vi2[over2inL2]=v2b
    u3b,v3b=cube.ca2cs_eqv2(us3,vs3,cos_o2inveclat3,cos_o2inveclon3,tan_o2inveclat3,tan_o2inveclon3);ui3[over2inL3]=u3b;vi3[over2inL3]=v3b
    u4b,v4b=cube.ca2cs_eqv2(us4,vs4,cos_o2inveclat4,cos_o2inveclon4,tan_o2inveclat4,tan_o2inveclon4);ui4[over2inL4]=u4b;vi4[over2inL4]=v4b
    u5b,v5b=cube.ca2cs_npv2(us5,vs5,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);ui5[over2inL5]=u5b;vi5[over2inL5]=v5b
    u6b,v6b=cube.ca2cs_spv2(us6,vs6,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);ui6[over2inL6]=u6b;vi6[over2inL6]=v6b

    cu1b,cv1b=cube.calc_covariant_eqv2(u1b,v1b,cos_o2inveclat1,cos_o2inveclon1,sin_o2inveclat1,sin_o2inveclon1);cut1[over2inL1]=cu1b;cvt1[over2inL1]=cv1b
    cu2b,cv2b=cube.calc_covariant_eqv2(u2b,v2b,cos_o2inveclat2,cos_o2inveclon2,sin_o2inveclat2,sin_o2inveclon2);cut2[over2inL2]=cu2b;cvt2[over2inL2]=cv2b
    cu3b,cv3b=cube.calc_covariant_eqv2(u3b,v3b,cos_o2inveclat3,cos_o2inveclon3,sin_o2inveclat3,sin_o2inveclon3);cut3[over2inL3]=cu3b;cvt3[over2inL3]=cv3b
    cu4b,cv4b=cube.calc_covariant_eqv2(u4b,v4b,cos_o2inveclat4,cos_o2inveclon4,sin_o2inveclat4,sin_o2inveclon4);cut4[over2inL4]=cu4b;cvt4[over2inL4]=cv4b
    cu5b,cv5b=cube.calc_covariant_npv2(u5b,v5b,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);cut5[over2inL5]=cu5b;cvt5[over2inL5]=cv5b
    cu6b,cv6b=cube.calc_covariant_spv2(u6b,v6b,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);cut6[over2inL6]=cu6b;cvt6[over2inL6]=cv6b
    
    # calc_tendency
    du1,dv1,dh1=integration.calc_tendencyv2(SDx1,SDy1,hiv1,cut1,cvt1,ui1,vi1,ht1,gi1,col1,dt)
    du2,dv2,dh2=integration.calc_tendencyv2(SDx2,SDy2,hiv2,cut2,cvt2,ui2,vi2,ht2,gi2,col2,dt)
    du3,dv3,dh3=integration.calc_tendencyv2(SDx3,SDy3,hiv3,cut3,cvt3,ui3,vi3,ht3,gi3,col3,dt)
    du4,dv4,dh4=integration.calc_tendencyv2(SDx4,SDy4,hiv4,cut4,cvt4,ui4,vi4,ht4,gi4,col4,dt)
    du5,dv5,dh5=integration.calc_tendencyv2(SDx5,SDy5,hiv5,cut5,cvt5,ui5,vi5,ht5,gi5,col5,dt)
    du6,dv6,dh6=integration.calc_tendencyv2(SDx6,SDy6,hiv6,cut6,cvt6,ui6,vi6,ht6,gi6,col6,dt)

    
    
    hi1+=dh1/3;hi3+=dh3/3
    hi4+=dh4/3;hi2+=dh2/3
    hi5+=dh5/3;hi6+=dh6/3
    cui1+=du1/3;cui2+=du2/3;cui3+=du3/3
    cui4+=du4/3;cui5+=du5/3;cui6+=du6/3
    cvi1+=dv1/3;cvi2+=dv2/3;cvi3+=dv3/3
    cvi4+=dv4/3;cvi5+=dv5/3;cvi6+=dv6/3

    
    #K4
    cut1=cu1c+du1;cut2=cu2c+du2;cut3=cu3c+du3
    cut4=cu4c+du4;cut5=cu5c+du5;cut6=cu6c+du6

    cvt1=cv1c+dv1;cvt2=cv2c+dv2;cvt3=cv3c+dv3
    cvt4=cv4c+dv4;cvt5=cv5c+dv5;cvt6=cv6c+dv6

    ht1=h1c+dh1;ht2=h2c+dh2;ht3=h3c+dh3
    ht4=h4c+dh4;ht5=h5c+dh5;ht6=h6c+dh6
    
    #input vals from planevals to overlapvector    
    hb1=ht1[in2bR1];ho[b2overlap1]=hb1.copy();hb2=ht2[in2bR2];ho[b2overlap2]=hb2.copy()
    hb3=ht3[in2bR3];ho[b2overlap3]=hb3.copy();hb4=ht4[in2bR4];ho[b2overlap4]=hb4.copy()
    hb5=ht5[in2bR5];ho[b2overlap5]=hb5.copy();hb6=ht6[in2bR6];ho[b2overlap6]=hb6.copy()

    
     
    # calc contravariant components 
    ui1,vi1=cube.co2contra_eqv2(cut1,cvt1,coslat1,coslon1,sinlat1,sinlon1)
    ui2,vi2=cube.co2contra_eqv2(cut2,cvt2,coslat2,coslon2,sinlat2,sinlon2)
    ui3,vi3=cube.co2contra_eqv2(cut3,cvt3,coslat3,coslon3,sinlat3,sinlon3)
    ui4,vi4=cube.co2contra_eqv2(cut4,cvt4,coslat4,coslon4,sinlat4,sinlon4)
    ui5,vi5=cube.co2contra_npv2(cut5,cvt5,coslat5,coslon5,sinlat5,sinlon5)
    ui6,vi6=cube.co2contra_spv2(cut6,cvt6,coslat6,coslon6,sinlat6,sinlon6)
    u1b=ui1[in2bR1];u2b=ui2[in2bR2];u3b=ui3[in2bR3];u4b=ui4[in2bR4];u5b=ui5[in2bR5];u6b=ui6[in2bR6]
    v1b=vi1[in2bR1];v2b=vi2[in2bR2];v3b=vi3[in2bR3];v4b=vi4[in2bR4];v5b=vi5[in2bR5];v6b=vi6[in2bR6]
      
    
    #transform  vector form contravariant to sphere
    us1,vs1=cube.cs2ca_eqv2(u1b,v1b,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
    us2,vs2=cube.cs2ca_eqv2(u2b,v2b,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
    us3,vs3=cube.cs2ca_eqv2(u3b,v3b,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
    us4,vs4=cube.cs2ca_eqv2(u4b,v4b,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
    us5,vs5=cube.cs2ca_npv2(u5b,v5b,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
    us6,vs6=cube.cs2ca_spv2(u6b,v6b,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)
    

    us[b2overlap1]=us1.copy();us[b2overlap2]=us2.copy();us[b2overlap3]=us3.copy();us[b2overlap4]=us4.copy();us[b2overlap5]=us5.copy();us[b2overlap6]=us6.copy()
    vs[b2overlap1]=vs1.copy();vs[b2overlap2]=vs2.copy();vs[b2overlap3]=vs3.copy();vs[b2overlap4]=vs4.copy();vs[b2overlap5]=vs5.copy();vs[b2overlap6]=vs6.copy()
    
    
    # overlaparea to difvector
    ht1[over2inL1]=ho[over2inR1];ht2[over2inL2]=ho[over2inR2];ht3[over2inL3]=ho[over2inR3]
    ht4[over2inL4]=ho[over2inR4];ht5[over2inL5]=ho[over2inR5];ht6[over2inL6]=ho[over2inR6]


    us1=us[over2inR1];us2=us[over2inR2];us3=us[over2inR3];us4=us[over2inR4];us5=us[over2inR5];us6=us[over2inR6]
    vs1=vs[over2inR1];vs2=vs[over2inR2];vs3=vs[over2inR3];vs4=vs[over2inR4];vs5=vs[over2inR5];vs6=vs[over2inR6]
    u1b,v1b=cube.ca2cs_eqv2(us1,vs1,cos_o2inveclat1,cos_o2inveclon1,tan_o2inveclat1,tan_o2inveclon1);ui1[over2inL1]=u1b;vi1[over2inL1]=v1b
    u2b,v2b=cube.ca2cs_eqv2(us2,vs2,cos_o2inveclat2,cos_o2inveclon2,tan_o2inveclat2,tan_o2inveclon2);ui2[over2inL2]=u2b;vi2[over2inL2]=v2b
    u3b,v3b=cube.ca2cs_eqv2(us3,vs3,cos_o2inveclat3,cos_o2inveclon3,tan_o2inveclat3,tan_o2inveclon3);ui3[over2inL3]=u3b;vi3[over2inL3]=v3b
    u4b,v4b=cube.ca2cs_eqv2(us4,vs4,cos_o2inveclat4,cos_o2inveclon4,tan_o2inveclat4,tan_o2inveclon4);ui4[over2inL4]=u4b;vi4[over2inL4]=v4b
    u5b,v5b=cube.ca2cs_npv2(us5,vs5,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);ui5[over2inL5]=u5b;vi5[over2inL5]=v5b
    u6b,v6b=cube.ca2cs_spv2(us6,vs6,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);ui6[over2inL6]=u6b;vi6[over2inL6]=v6b

    cu1b,cv1b=cube.calc_covariant_eqv2(u1b,v1b,cos_o2inveclat1,cos_o2inveclon1,sin_o2inveclat1,sin_o2inveclon1);cut1[over2inL1]=cu1b;cvt1[over2inL1]=cv1b
    cu2b,cv2b=cube.calc_covariant_eqv2(u2b,v2b,cos_o2inveclat2,cos_o2inveclon2,sin_o2inveclat2,sin_o2inveclon2);cut2[over2inL2]=cu2b;cvt2[over2inL2]=cv2b
    cu3b,cv3b=cube.calc_covariant_eqv2(u3b,v3b,cos_o2inveclat3,cos_o2inveclon3,sin_o2inveclat3,sin_o2inveclon3);cut3[over2inL3]=cu3b;cvt3[over2inL3]=cv3b
    cu4b,cv4b=cube.calc_covariant_eqv2(u4b,v4b,cos_o2inveclat4,cos_o2inveclon4,sin_o2inveclat4,sin_o2inveclon4);cut4[over2inL4]=cu4b;cvt4[over2inL4]=cv4b
    cu5b,cv5b=cube.calc_covariant_npv2(u5b,v5b,cos_o2inveclat5,cos_o2inveclon5,sin_o2inveclat5,sin_o2inveclon5);cut5[over2inL5]=cu5b;cvt5[over2inL5]=cv5b
    cu6b,cv6b=cube.calc_covariant_spv2(u6b,v6b,cos_o2inveclat6,cos_o2inveclon6,sin_o2inveclat6,sin_o2inveclon6);cut6[over2inL6]=cu6b;cvt6[over2inL6]=cv6b
    
    # calc_tendency
    du1,dv1,dh1=integration.calc_tendencyv2(SDx1,SDy1,hiv1,cut1,cvt1,ui1,vi1,ht1,gi1,col1,dt)
    du2,dv2,dh2=integration.calc_tendencyv2(SDx2,SDy2,hiv2,cut2,cvt2,ui2,vi2,ht2,gi2,col2,dt)
    du3,dv3,dh3=integration.calc_tendencyv2(SDx3,SDy3,hiv3,cut3,cvt3,ui3,vi3,ht3,gi3,col3,dt)
    du4,dv4,dh4=integration.calc_tendencyv2(SDx4,SDy4,hiv4,cut4,cvt4,ui4,vi4,ht4,gi4,col4,dt)
    du5,dv5,dh5=integration.calc_tendencyv2(SDx5,SDy5,hiv5,cut5,cvt5,ui5,vi5,ht5,gi5,col5,dt)
    du6,dv6,dh6=integration.calc_tendencyv2(SDx6,SDy6,hiv6,cut6,cvt6,ui6,vi6,ht6,gi6,col6,dt)

    
    hi1+=dh1/6;hi3+=dh3/6
    hi4+=dh4/6;hi2+=dh2/6
    hi5+=dh5/6;hi6+=dh6/6
    cui1+=du1/6;cui2+=du2/6;cui3+=du3/6
    cui4+=du4/6;cui5+=du5/6;cui6+=du6/6
    cvi1+=dv1/6;cvi2+=dv2/6;cvi3+=dv3/6
    cvi4+=dv4/6;cvi5+=dv5/6;cvi6+=dv6/6
    
    #input vals in palne to overlap

    hb1=hi1[in2bR1];ho[b2overlap1]=hb1;hb2=hi2[in2bR2];ho[b2overlap2]=hb2
    hb3=hi3[in2bR3];ho[b2overlap3]=hb3;hb4=hi4[in2bR4];ho[b2overlap4]=hb4
    hb5=hi5[in2bR5];ho[b2overlap5]=hb5;hb6=hi6[in2bR6];ho[b2overlap6]=hb6

    u1b,v1b=cube.co2contra_eqv2(cui1[in2bR1],cvi1[in2bR1],cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
    u2b,v2b=cube.co2contra_eqv2(cui2[in2bR2],cvi2[in2bR2],cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
    u3b,v3b=cube.co2contra_eqv2(cui3[in2bR3],cvi3[in2bR3],cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
    u4b,v4b=cube.co2contra_eqv2(cui4[in2bR4],cvi4[in2bR4],cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
    u5b,v5b=cube.co2contra_npv2(cui5[in2bR5],cvi5[in2bR5],cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
    u6b,v6b=cube.co2contra_spv2(cui6[in2bR6],cvi6[in2bR6],cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)    

    us1,vs1=cube.cs2ca_eqv2(u1b,v1b,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
    us2,vs2=cube.cs2ca_eqv2(u2b,v2b,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
    us3,vs3=cube.cs2ca_eqv2(u3b,v3b,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
    us4,vs4=cube.cs2ca_eqv2(u4b,v4b,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
    us5,vs5=cube.cs2ca_npv2(u5b,v5b,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
    us6,vs6=cube.cs2ca_spv2(u6b,v6b,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)    

    us[b2overlap1]=us1;us[b2overlap2]=us2;us[b2overlap3]=us3
    us[b2overlap4]=us4;us[b2overlap5]=us5;us[b2overlap6]=us6
    vs[b2overlap1]=vs1;vs[b2overlap2]=vs2;vs[b2overlap3]=vs3;vs[b2overlap4]=vs4;vs[b2overlap5]=vs5;vs[b2overlap6]=vs6

    if i==srec-1:
        h[cidim1]=hi1[c2innerL1];h[cidim2]=hi2[c2innerL2];h[cidim3]=hi3[c2innerL3]
        h[cidim4]=hi4[c2innerL4];h[cidim5]=hi5[c2innerL5];h[cidim6]=hi6[c2innerL6]
        log_h[count,:]+=h
        acu[cidim1]=cui1[c2innerL1];acv[cidim1]=cvi1[c2innerL1]
        acu[cidim2]=cui2[c2innerL2];acv[cidim2]=cvi2[c2innerL2]
        acu[cidim3]=cui3[c2innerL3];acv[cidim3]=cvi3[c2innerL3]
        acu[cidim4]=cui4[c2innerL4];acv[cidim4]=cvi4[c2innerL4]
        acu[cidim5]=cui5[c2innerL5];acv[cidim5]=cvi5[c2innerL5]
        acu[cidim6]=cui6[c2innerL6];acv[cidim6]=cvi6[c2innerL6]
        log_cu[count,:]=acu
        log_cv[count,:]=acv

        srec+=srec0
        count+=1

    


endtime=time.time()
print('the integration is finished in {:.1f}s'.format(endtime-start))

#np.save('case{}_alpha{:d}_log_h_{}N{}n{}deg{}K{}dt{}gamma{}hep{}'.format(case,int(alpha/pi*180),node,N,n,deg,k,dt//60,gamma,hep),log_h)
#np.save('case{}_alpha{:d}_log_cu_{}N{}n{}deg{}K{}dt{}gamma{}hep{}'.format(case,int(alpha/pi*180),node,N,n,deg,k,dt//60,gamma,hep),log_cu)
#np.save('case{}_alpha{:d}_log_cv_{}N{}n{}deg{}K{}dt{}gamma{}hep{}'.format(case,int(alpha/pi*180),node,N,n,deg,k,dt//60,gamma,hep),log_cv)

triang = tri.Triangulation(lon/pi*180, lat/pi*180)
fig=plt.figure(figsize=(12,6))
ax=fig.add_subplot(111)
if case==2:
    lev=11#np.linspace(-1e-3,1e-3,11,endpoint=True)
    tcf=ax.tricontourf(triang,h-log_h[0],levels=lev,cmap='coolwarm')
    tcfb=ax.tricontour(triang,h-log_h[0],levels=lev,colors='k')
    plt.colorbar(tcf,ax=ax)
    #plt.savefig('case2_alpha{:d}_error_h_{}N{}n{}deg{}K{}dt{}gamma{}hep{}.png'.format(int(alpha/pi*180),node,N,n,deg,k,dt//60,gamma,hep))
elif case==6:
    tcf=ax.tricontourf(triang,h,levels=np.linspace(8000,10800,15,endpoint=True),cmap='jet')
    tcfb=ax.tricontour(triang,h,levels=np.linspace(8000,10800,15,endpoint=True),colors='k')
    plt.colorbar(tcf,ax=ax)
    plt.savefig('case6_{}N{}n{}deg{}K{}dt{}gamma{}hep{}.png'.format(node,N,n,deg,k,dt//60,gamma,hep))

plt.show()


# calc error norm
if case==2:
    ht=log_h[0]
    h0=wgt@abs(log_h[0,:])
    h1=wgt@abs(log_h[0,:]-h)
    l1=h1/h0
    h0h0=(wgt@(ht**2))**(1/2)
    dh=(wgt@((h-ht)**2))**(1/2)
    l2=dh/h0h0
    linf=np.amax(abs(log_h[0,:]-h))/np.amax(abs(log_h[0]))

    print('l1={:e}'.format(l1))
    print('l2={:e}'.format(l2))
    print('linf={:e}'.format(linf))
    hozon=((wgt@h)/(wgt@ht)-1)
    print('conservation errro={:.1e}'.format(hozon))

