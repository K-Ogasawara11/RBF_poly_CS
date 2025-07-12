import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from numpy import pi,sin,cos
import scipy
import initial_condition as init
import matplotlib.tri as tri
from tqdm import tqdm
import md_cube as cube
import md_dmatrix as dmat
import os
import time
import md_adv as mdadv
from params import N,node,n,hep,deg,d,k,gamma,bell,day,step_in_day,rec_day,alpha,lon0,lat0

ep=hep

#the coefficient of hyperviscosity for transport sheme (Gumderman et al. 2020)
delta=-(1600/N)**k*gamma;print('delta{:.1e}'.format(delta))
stability_domain='False'
a=1/np.sqrt(3)

# parameter for transport scheme (don't change)
projection='edist'
dt=1

#定数の設定
a=1/np.sqrt(3)

# load nodes

if node=='md':
    lon,lat,x,y,z,wgt=cube.knn_md(N)
elif node=='icos':
    lon,lat,x,y,z=cube.knn_icos(N)
    wgt=np.loadtxt('knn_icos_nomodify_weight_N{}.txt'.format(N))



print('np.amax(lat)',np.amax(lat))
#風（極座標）
alpha=pi/4
lat0=0
lon0=-pi/2
u0=2*pi/day
#u=u0*np.cos(lat)
u=u0*(np.cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
#v=np.zeros(N)
v=-u0*sin(lon)*sin(alpha)

#トレーサーの計算。


if bell=='cos':
    h =init.cosbell2(lon, lat,lon0,lat0,u0,alpha,0.0,1,1/3)
elif bell=='gauss': 
    h=init.gaussbell2(lon, lat,lon0,lat0,u0,alpha,0.0,1,1/3)
elif bell=='gund':
    h=init.gaussbell3(x,lon0,lon,lat)


idim1,idim2,idim3,idim4,idim5,idim6=cube.partition_xyz(x,y,z,lat,lon,N)
p1_vec2all,p1_all2vec=mdadv.make_transform_matrix(idim1,N)
p2_vec2all,p2_all2vec=mdadv.make_transform_matrix(idim2,N)
p3_vec2all,p3_all2vec=mdadv.make_transform_matrix(idim3,N)
p4_vec2all,p4_all2vec=mdadv.make_transform_matrix(idim4,N)
p5_vec2all,p5_all2vec=mdadv.make_transform_matrix(idim5,N)
p6_vec2all,p6_all2vec=mdadv.make_transform_matrix(idim6,N)

x1=p1_all2vec@x;y1=p1_all2vec@y;z1=p1_all2vec@z
x2=p2_all2vec@x;y2=p2_all2vec@y;z2=p2_all2vec@z
x3=p3_all2vec@x;y3=p3_all2vec@y;z3=p3_all2vec@z
x4=p4_all2vec@x;y4=p4_all2vec@y;z4=p4_all2vec@z
x5=p5_all2vec@x;y5=p5_all2vec@y;z5=p5_all2vec@z
x6=p6_all2vec@x;y6=p6_all2vec@y;z6=p6_all2vec@z


lon1=p1_all2vec@lon;lat1=p1_all2vec@lat;g1=cos(lat1)**3*cos(lon1)**3/a**2
lon2=p2_all2vec@lon;lat2=p2_all2vec@lat
lon2=lon2-pi/2
g2=cos(lat2)**3*cos(lon2)**3/a**2
lon3=p3_all2vec@lon;lat3=p3_all2vec@lat
for i in range(lon3.size):
    if lon3[i]<0:
        lon3[i]=lon3[i]+pi*2
lon3-=pi
g3=cos(lat3)**3*cos(lon3)**3/a**2
lon4=p4_all2vec@lon;lat4=p4_all2vec@lat
lon4+=pi/2
g4=cos(lat4)**3*cos(lon4)**3/a**2
lon5=p5_all2vec@lon;lat5=p5_all2vec@lat;g5=sin(lat5)**3/a**2
lon6=p6_all2vec@lon;lat6=p6_all2vec@lat;g6=sin(lat6)**3/a**2



gvec=np.zeros(N)
gvec+=p1_vec2all@g1
gvec+=p2_vec2all@g2
gvec+=p3_vec2all@g3
gvec+=p4_vec2all@g4
gvec+=p5_vec2all@g5
gvec+=p6_vec2all@g6

#calc gh
gh=h.copy()

if os.path.isfile('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx1=np.load('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy1=np.load('../Dmatrix/Dy1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index1=np.load('../Dmatrix/index1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx1,Dy1,r,index1=mdadv.calc_Dmatrix_phs(x,y,z,x1,y1,z1,idim1,1,n,deg,projection)
    np.save('../Dmatrix/Dx1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx1)
    np.save('../Dmatrix/Dy1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy1)
    np.save('../Dmatrix/index1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index1)
Dx1,Dy1=mdadv.metric_eq_flux(Dx1,Dy1,index1,idim1,u,v,lat,lon,n,1)

#P2

if os.path.isfile('../Dmatrix/Dx2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx2=np.load('../Dmatrix/Dx2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy2=np.load('../Dmatrix/Dy2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index2=np.load('../Dmatrix/index2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx2,Dy2,r,index2=mdadv.calc_Dmatrix_phs(x,y,z,x2,y2,z2,idim2,2,n,deg,projection)
    np.save('../Dmatrix/Dx2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx2)
    np.save('../Dmatrix/Dy2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy2)
    np.save('../Dmatrix/index2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index2)

Dx2,Dy2=mdadv.metric_eq_flux(Dx2,Dy2,index2,idim2,u,v,lat,lon,n,2)



#P3

if os.path.isfile('../Dmatrix/Dx3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx3=np.load('../Dmatrix/Dx3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy3=np.load('../Dmatrix/Dy3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index3=np.load('../Dmatrix/index3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx3,Dy3,r,index3=mdadv.calc_Dmatrix_phs(x,y,z,x3,y3,z3,idim3,3,n,deg,projection)
    np.save('../Dmatrix/Dx3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx3)
    np.save('../Dmatrix/Dy3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy3)
    np.save('../Dmatrix/index3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index3)
Dx3,Dy3=mdadv.metric_eq_flux(Dx3,Dy3,index3,idim3,u,v,lat,lon,n,3)
#SDx3,SDy3=dmat.convert_sparse(index,Dx,Dy,idim3,N,n)

#P4

if os.path.isfile('../Dmatrix/Dx4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx4=np.load('../Dmatrix/Dx4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy4=np.load('../Dmatrix/Dy4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index4=np.load('../Dmatrix/index4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx4,Dy4,r,index4=mdadv.calc_Dmatrix_phs(x,y,z,x4,y4,z4,idim4,4,n,deg,projection)
    np.save('../Dmatrix/Dx4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx4)
    np.save('../Dmatrix/Dy4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy4)
    np.save('../Dmatrix/index4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index4)
Dx4,Dy4=mdadv.metric_eq_flux(Dx4,Dy4,index4,idim4,u,v,lat,lon,n,4)

#P5
if os.path.isfile('../Dmatrix/Dx5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx5=np.load('../Dmatrix/Dx5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy5=np.load('../Dmatrix/Dy5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index5=np.load('../Dmatrix/index5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx5,Dy5,r,index5=mdadv.calc_Dmatrix_phs(x,y,z,x5,y5,z5,idim5,5,n,deg,projection)
    np.save('../Dmatrix/Dx5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx5)
    np.save('../Dmatrix/Dy5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy5)
    np.save('../Dmatrix/index5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index5)
Dx5,Dy5=mdadv.metric_np_flux(Dx5,Dy5,index5,idim5,u,v,lat,lon,n)

#P6
if os.path.isfile('../Dmatrix/Dx6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx6=np.load('../Dmatrix/Dx6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy6=np.load('../Dmatrix/Dy6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index6=np.load('../Dmatrix/index6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx6,Dy6,r,index6=mdadv.calc_Dmatrix_phs(x,y,z,x6,y6,z6,idim6,6,n,deg,projection)
    np.save('../Dmatrix/Dx6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx6)
    np.save('../Dmatrix/Dy6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy6)
    np.save('../Dmatrix/index6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index6)
Dx6,Dy6=mdadv.metric_sp_flux(Dx6,Dy6,index6,idim6,u,v,lat,lon,n)



if os.path.isfile('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))==True:
    H=np.load('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
    indexb=np.load('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
else:
    indexb=mdadv.calc_index_pc(N,n,x,y,z)
    R=mdadv.calc_R_pc(N,n,x,y,z,indexb)
    H=mdadv.calc_hiv_FH2011(d,k,hep,R,N,n)
    np.save('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),indexb)
#row,col,vals=mdadv.make_row_col_pc(indexb,H*delta,n,N)
#hiv=scipy.sparse.csr_matrix((vals,(row,col)),shape = (N,N))


print(N,n,N*n,idim1.shape)
D=np.zeros((N,N))


count=0
for i in range(idim1.size):
    indi1=int(idim1[i])
    for j in range(n):
        indj1=int(index1[i,j])
        D[indi1,indj1]=Dx1[i,j]+Dy1[i,j]
        count+=1
count=0
for i in range(idim2.size):
    indi1=int(idim2[i])
    for j in range(n):
        indj1=int(index2[i,j])
        D[indi1,indj1]=Dx2[i,j]+Dy2[i,j]
        count+=1
count=0
for i in range(idim3.size):
    indi1=int(idim3[i])
    for j in range(n):
        indj1=int(index3[i,j])
        D[indi1,indj1]=Dx3[i,j]+Dy3[i,j]
        count+=1
count=0
for i in range(idim4.size):
    indi1=int(idim4[i])
    for j in range(n):
        indj1=int(index4[i,j])
        D[indi1,indj1]=Dx4[i,j]+Dy4[i,j]
        count+=1

count=0
for i in range(idim5.size):
    indi1=int(idim5[i])
    for j in range(n):
        indj1=int(index5[i,j])
        D[indi1,indj1]=Dx5[i,j]+Dy5[i,j]
        count+=1

count=0
for i in range(idim6.size):
    indi1=int(idim6[i])
    for j in range(n):
        indj1=int(index6[i,j])
        D[indi1,indj1]=Dx6[i,j]+Dy6[i,j]
        count+=1


count=0
for i in range(N):
    indi1=int(indexb[i,0])
    for j in range(n):
        indj1=int(indexb[i,j])
        D[indi1,indj1]+=H[i,j]*delta
        count+=1


eigv=np.linalg.eigvals(D)
if gamma != 0:
    np.save('fluxform_eigvals_{}N{}n{}deg{}_with_hivK{}_gamm{:.1e}_hep{}'.format(node,N,n,deg,k,gamma,hep),eigv)
else:
    np.save('fluxform_eigvals_{}N{}n{}deg{}_no_hiv_hep{}'.format(node,N,n,deg,hep),eigv)


reigv=np.real(eigv)
ieigv=np.imag(eigv)
plt.scatter(reigv,ieigv)
maxreig=np.amax(reigv[reigv>0])
print(maxreig)

#points=3000
#x2=np.linspace(-3,0.5,points,endpoint=True)
#y2=np.linspace(-3,3,points,endpoint=True)
#X,Y=np.meshgrid(x2,y2)
#z=X+Y*1j
#rz=np.real(np.absolute(1+z+(z**2)/2+(z**3)/6+(z**4)/24))
#for i in range(points):
#    for j in range(points):
#        if rz[i,j]>1:
#            rz[i,j]=0
#        elif rz[i,j]<0.999:
#            rz[i,j]=0
#plt.contour(X,Y,rz,colors='k',linewidths=1)


plt.xlabel('real')
plt.ylabel('image')
plt.title('eigenvalues dt={}min'.format(dt//60))

print('max real eigvals',np.amax(reigv))

#plt.savefig('eigenvalue_{}_N{}n{}dt{}_withRK4_absolute_stability_with_hiv.jpeg'.format(node,N,n,dt//60))
#plt.savefig('eigenvalue_{}_N{}n{}dt{}_with_hivK{}_hep{}gamma{:.1e}.jpeg'.format(node,N,n,dt//60,k,hep,gamma))
plt.show()