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
import md_integration as integration
import os
import sys
import md_adv as mdadv
from params import N,node,n,hep,deg,d,k,gamma,bell,day,step_in_day,rec_day,alpha,lon0,lat0

#the coefficient of hyperviscosity for transport sheme (Gumderman et al. 2020)
delta=-(1600/N)**k*gamma;print('delta{:.1e}'.format(delta))


# parameter for transport scheme (don't change)
projection='edist'
rot=day/12
cday=12*rot
steps=int(step_in_day*cday)
dt=1/(step_in_day);print('dt~{}min'.format(24*60//step_in_day))
nrec=int(cday/rec_day);print('rec no. {}'.format(nrec))
srec=step_in_day*rec_day
srec0=step_in_day*rec_day
log_h=np.zeros((nrec+1,N))

#params
a=1/np.sqrt(3)

# load nodes

if node=='md':
    lon,lat,x,y,z,wgt=cube.knn_md(N)
elif node=='icos':
    lon,lat,x,y,z=cube.knn_icos(N)
    wgt=np.loadtxt('knn_icos_nomodify_weight_N{}.txt'.format(N))

elif node=='bauer':
    lon,lat,x,y,z=cube.bauer(N)
    wgt=np.loadtxt('../nodes/bauer_weight_N{}.txt'.format(N)) 


# initialization for transport sheme
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

    # set or calc Dmatrix

#P1
if os.path.isfile('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index1=np.load('../Dmatrix/index1_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index1=mdadv.calc_Dmatrix_phs(x,y,z,x1,y1,z1,idim1,1,n,deg,projection)
    np.save('../Dmatrix/Dx1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index1_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index1)
print(idim1.shape,Dx.shape)
SDx1,SDy1=mdadv.convert_sparse(index1,Dx,Dy,idim1,N,n)


#P2

if os.path.isfile('../Dmatrix/Dx2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index2=np.load('../Dmatrix/index2_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index2=mdadv.calc_Dmatrix_phs(x,y,z,x2,y2,z2,idim2,2,n,deg,projection)
    np.save('../Dmatrix/Dx2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index2_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index2)

SDx2,SDy2=mdadv.convert_sparse(index2,Dx,Dy,idim2,N,n)


#P3

if os.path.isfile('../Dmatrix/Dx3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index3=np.load('../Dmatrix/index3_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index3=mdadv.calc_Dmatrix_phs(x,y,z,x3,y3,z3,idim3,3,n,deg,projection)
    np.save('../Dmatrix/Dx3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index3_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index3)

SDx3,SDy3=mdadv.convert_sparse(index3,Dx,Dy,idim3,N,n)

#P4

if os.path.isfile('../Dmatrix/Dx4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index4=np.load('../Dmatrix/index4_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index4=mdadv.calc_Dmatrix_phs(x,y,z,x4,y4,z4,idim4,4,n,deg,projection)
    np.save('../Dmatrix/Dx4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index4_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index4)
SDx4,SDy4=mdadv.convert_sparse(index4,Dx,Dy,idim4,N,n)

#P5

if os.path.isfile('../Dmatrix/Dx5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index5=np.load('../Dmatrix/index5_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index5=mdadv.calc_Dmatrix_phs(x,y,z,x5,y5,z5,idim5,5,n,deg,projection)
    np.save('../Dmatrix/Dx5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index5_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index5)
SDx5,SDy5=mdadv.convert_sparse(index5,Dx,Dy,idim5,N,n)

#P6
if os.path.isfile('../Dmatrix/Dx6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))==True:
    Dx=np.load('../Dmatrix/Dx6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    Dy=np.load('../Dmatrix/Dy6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
    index6=np.load('../Dmatrix/index6_adv_{}N{}n{}deg{}.npy'.format(node,N,n,deg))
else:
    Dx,Dy,r,index6=mdadv.calc_Dmatrix_phs(x,y,z,x6,y6,z6,idim6,6,n,deg,projection)
    np.save('../Dmatrix/Dx6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dx)
    np.save('../Dmatrix/Dy6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),Dy)
    np.save('../Dmatrix/index6_adv_{}N{}n{}deg{}'.format(node,N,n,deg),index6)

SDx6,SDy6=mdadv.convert_sparse(index6,Dx,Dy,idim6,N,n)



if os.path.isfile('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))==True:
    H=np.load('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
    indexb=np.load('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}.npy'.format(node,N,n,hep,k))
else:
    indexb=mdadv.calc_index_pc(N,n,x,y,z)
    R=mdadv.calc_R_pc(N,n,x,y,z,indexb)
    H=mdadv.calc_hiv_FH2011(d,k,hep,R,N,n)
    np.save('../Dmatrix/HFL_adv_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),H)
    np.save('../Dmatrix/HFL_adv_index_{}N{}n{}ep{:.1f}k{}'.format(node,N,n,hep,k),indexb)
row,col,vals=mdadv.make_row_col_pc(indexb,H*delta,n,N)
hiv=scipy.sparse.csr_matrix((vals,(row,col)),shape = (N,N))


# calc initial condition
u0=2*pi/12
u=u0*(np.cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
v=-u0*sin(lon)*sin(alpha)

if bell=='cos':
    h =init.cosbell2(lon, lat,lon0,lat0,u0,alpha,0.0,1,1/3)
elif bell=='gauss': 
    h=init.gaussbell2(lon, lat,lon0,lat0,u0,alpha,0.0,1,1/3)
elif bell=='gund':
    h=init.gaussbell3(x,lon0,lon,lat)

gh=h.copy()
log_h[0,:]=h.copy()

u1,v1=mdadv.ca2cs_eq(u[idim1],v[idim1],lat[idim1],lon[idim1],1)
u2,v2=mdadv.ca2cs_eq(u[idim2],v[idim2],lat[idim2],lon[idim2],2)
u3,v3=mdadv.ca2cs_eq(u[idim3],v[idim3],lat[idim3],lon[idim3],3)
u4,v4=mdadv.ca2cs_eq(u[idim4],v[idim4],lat[idim4],lon[idim4],4)
u5,v5=mdadv.ca2cs_np(u[idim5],v[idim5],lat[idim5],lon[idim5])
u6,v6=mdadv.ca2cs_sp(u[idim6],v[idim6],lat[idim6],lon[idim6])



#RK4
dh=np.zeros(N)
count=1
triang = tri.Triangulation(lon, lat)
fig = plt.figure(figsize=(7,7))
for i in tqdm(range(steps)):
    hc=gh.copy()
    #K1
    dhb=integration.calc_dh_adv(SDx1,SDy1,u1,v1,p1_vec2all,hc)
    dh=-dhb
    dhb=integration.calc_dh_adv(SDx2,SDy2,u2,v2,p2_vec2all,hc)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx3,SDy3,u3,v3,p3_vec2all,hc)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx4,SDy4,u4,v4,p4_vec2all,hc)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx5,SDy5,u5,v5,p5_vec2all,hc)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx6,SDy6,u6,v6,p6_vec2all,hc)
    dh-=dhb
    dh+=(hiv@hc)
    dh=dh*dt
    gh+=dh/6

    #K2
    ht=(hc+dh/2)
    dhb=integration.calc_dh_adv(SDx1,SDy1,u1,v1,p1_vec2all,ht)
    dh=-dhb
    dhb=integration.calc_dh_adv(SDx2,SDy2,u2,v2,p2_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx3,SDy3,u3,v3,p3_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx4,SDy4,u4,v4,p4_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx5,SDy5,u5,v5,p5_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx6,SDy6,u6,v6,p6_vec2all,ht)
    dh-=dhb
    dh+=(hiv@ht)
    dh=dh*dt
    gh+=dh/3
    

    #K3
    ht=(hc+dh/2)
    
    dhb=integration.calc_dh_adv(SDx1,SDy1,u1,v1,p1_vec2all,ht)
    dh=-dhb
    dhb=integration.calc_dh_adv(SDx2,SDy2,u2,v2,p2_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx3,SDy3,u3,v3,p3_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx4,SDy4,u4,v4,p4_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx5,SDy5,u5,v5,p5_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx6,SDy6,u6,v6,p6_vec2all,ht)
    dh-=dhb
    dh+=(hiv@ht)
    dh=dh*dt
    gh+=dh/3

    #K4
    ht=(hc+dh)
    
    dhb=integration.calc_dh_adv(SDx1,SDy1,u1,v1,p1_vec2all,ht)
    dh=-dhb
    dhb=integration.calc_dh_adv(SDx2,SDy2,u2,v2,p2_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx3,SDy3,u3,v3,p3_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx4,SDy4,u4,v4,p4_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx5,SDy5,u5,v5,p5_vec2all,ht)
    dh-=dhb
    dhb=integration.calc_dh_adv(SDx6,SDy6,u6,v6,p6_vec2all,ht)
    dh-=dhb
    dh+=(hiv@ht)
    dh=dh*dt
    gh+=dh/6

    if i==srec-1:
        log_h[count,:]+=gh
        srec+=srec0
        count+=1
        if np.amax(gh)>2e3:
            print('overflow')	
            sys.exit()

        
#np.savetxt('log_h_adv_days_{}bell_{}N{}n{}deg{}dt{:.1e}alpha{}_v2.txt'.format(bell,node,N,n,deg,dt,alpha/pi*180),log_h)
np.savetxt('{}rev_daycal_log_h_adv_{}bell_{}N{}n{}deg{}dt{}k{}gamma{}ep{}alpha45.txt'.format(rot,bell,node,N,n,deg,dt,k,gamma,hep),log_h)
triang = tri.Triangulation(lon, lat)
#fig = plt.figure(figsize=(7,7))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

tcf1 = ax1.tricontourf(triang,gh,levels=10,cmap=cm.bwr)
plt.colorbar(tcf1,ax=ax1)
tcf = ax2.tricontourf(triang,gh-log_h[0],levels=21,cmap=cm.jet)
#tcfb = ax2.tricontour(triang,log_h[0],levels=4,colors='k')
plt.colorbar(tcf,ax=ax2)
plt.show()


# calc error norm 
ht=log_h[0]
h0=wgt@abs(ht)
h1=wgt@abs(ht-gh);print(h1)
l1=h1/h0
h0h0=(wgt@(ht**2))**(1/2)
dh=(wgt@((gh-ht)**2))**(1/2)
l2=dh/h0h0
linf=np.amax(abs(log_h[0,:]-gh))/np.amax(abs(log_h[0]))
print('l1={:e}'.format(l1))
print('l2={:e}'.format(l2))
print('linf={:e}'.format(linf))
print('max',np.amax(gh))
print('min',np.amin(gh))



