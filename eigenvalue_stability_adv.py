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
import md_integration as integration
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

# load nodes

if node=='md':
    lon,lat,x,y,z,wgt=cube.knn_md(N)
elif node=='icos':
    lon,lat,x,y,z=cube.knn_icos(N)
    wgt=np.loadtxt('knn_icos_nomodify_weight_N{}.txt'.format(N))





# initialization for the cubed sphre projection
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



#calc initial condition
u0=2*pi/12
u=u0*(np.cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
v=-u0*sin(lon)*sin(alpha)

u1,v1=mdadv.ca2cs_eq(u[idim1],v[idim1],lat[idim1],lon[idim1],1)
u2,v2=mdadv.ca2cs_eq(u[idim2],v[idim2],lat[idim2],lon[idim2],2)
u3,v3=mdadv.ca2cs_eq(u[idim3],v[idim3],lat[idim3],lon[idim3],3)
u4,v4=mdadv.ca2cs_eq(u[idim4],v[idim4],lat[idim4],lon[idim4],4)
u5,v5=mdadv.ca2cs_np(u[idim5],v[idim5],lat[idim5],lon[idim5])
u6,v6=mdadv.ca2cs_sp(u[idim6],v[idim6],lat[idim6],lon[idim6])



#P1
if os.path.isfile('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx1=np.load('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy1=np.load('../Dmatrix/Dy1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index1=np.load('../Dmatrix/index1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx1,Dy1,r,index1=mdadv.calc_Dmatrix_phs(x,y,z,x1,y1,z1,idim1,1,n,deg,projection)
    np.save('../Dmatrix/Dx1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx1)
    np.save('../Dmatrix/Dy1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy1)
    np.save('../Dmatrix/index1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index1)


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




if os.path.isfile('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))==True:
    H=np.load('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
    indexb=np.load('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
else:
    indexb=mdadv.calc_index_pc(N,n,x,y,z)
    R=mdadv.calc_R_pc(N,n,x,y,z,indexb)
    H=mdadv.calc_hiv_FH2011(d,k,hep,R,N,n)
    np.save('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),indexb)





D=np.zeros((N,N))

for i in range(idim1.size):
    indi1=int(idim1[i])
    for j in range(n):
        indj1=int(index1[i,j])
        D[indi1,indj1]=u1[i]*Dx1[i,j]+v1[i]*Dy1[i,j]
        

for i in range(idim2.size):
    indi1=int(idim2[i])
    for j in range(n):
        indj1=int(index2[i,j])
        D[indi1,indj1]=u2[i]*Dx2[i,j]+v2[i]*Dy2[i,j]
        

for i in range(idim3.size):
    indi1=int(idim3[i])
    for j in range(n):
        indj1=int(index3[i,j])
        D[indi1,indj1]=u3[i]*Dx3[i,j]+v3[i]*Dy3[i,j]
        

for i in range(idim4.size):
    indi1=int(idim4[i])
    for j in range(n):
        indj1=int(index4[i,j])
        D[indi1,indj1]=u4[i]*Dx4[i,j]+v4[i]*Dy4[i,j]

for i in range(idim5.size):
    indi1=int(idim5[i])
    for j in range(n):
        indj1=int(index5[i,j])
        D[indi1,indj1]=u5[i]*Dx5[i,j]+v5[i]*Dy5[i,j]
        
for i in range(idim6.size):
    indi1=int(idim6[i])
    for j in range(n):
        indj1=int(index6[i,j])
        D[indi1,indj1]=u6[i]*Dx6[i,j]+v6[i]*Dy6[i,j]
        

for i in range(N):
    indi1=int(indexb[i,0])
    for j in range(n):
        indj1=int(indexb[i,j])
        D[indi1,indj1]+=H[i,j]*delta
   


print('start_calc_eigenvalue')
start=time.time()
eigv=np.linalg.eigvals(D)
if gamma != 0:
    np.save('advectionform_eigvals_{}N{}n{}deg{}_with_hivK{}_gamm{:.1e}_hep{}'.format(node,N,n,deg,k,gamma,hep),eigv)
else:
    np.save('advectionform_eigvals_{}N{}n{}deg{}_no_hiv_hep{}'.format(node,N,n,deg,hep),eigv)

endtime=time.time()
print('{:.1f}s'.format(endtime-start))
reigv=np.real(eigv)
ieigv=np.imag(eigv)
plt.scatter(reigv*dt,ieigv*dt)
#maxreig=np.amax(reigv[reigv>0])
#print(maxreig)
if stability_domain=='ture':

    points=3000
    x2=np.linspace(-3,0.5,points,endpoint=True)
    y2=np.linspace(-3,3,points,endpoint=True)
    X,Y=np.meshgrid(x2,y2)
    z=X+Y*1j
    rz=np.real(np.absolute(1+z+(z**2)/2+(z**3)/6+(z**4)/24))
    for i in range(points):
        for j in range(points):
            if rz[i,j]>1:
                rz[i,j]=0
            elif rz[i,j]<0.999:
                rz[i,j]=0
    plt.contour(X,Y,rz,colors='k',linewidths=1)

print(np.amin(abs(reigv)))
plt.xlabel('real')
plt.ylabel('image')
plt.title('eigenvalues dt={}min'.format(dt//60))

#print('max real eigvals',np.amax(reigv))

#plt.savefig('eigenvalue_{}_N{}n{}dt{}_withRK4_absolute_stability_with_hiv.jpeg'.format(node,N,n,dt//60))
#plt.savefig('eigenvalue_{}_N{}n{}dt{}.jpeg'.format(node,N,n,dt//60))
plt.show()