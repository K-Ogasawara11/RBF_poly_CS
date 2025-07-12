import numpy as np
import md_cube as cube
import md_ordering as order
from params import N,n,dt,steps,node,projection

print('N{},n{},dt{},steps{}'.format(N,n,dt,steps))
print('{} node'.format(node))
print('{} projection'.format(projection))

# load nodes

if node=='md':    
    lon,lat,x,y,z,wgt=cube.knn_md(N)

elif node=='icos_nomod':
    lon,lat,x,y,z=cube.knn_icos_nomod(N)
    wgt=np.loadtxt('knn_icos_nomodify_weight{}.txt'.format(N))


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

#output cindex
np.savetxt('cindex1.txt',cindex1,fmt='%d');np.savetxt('cindex2.txt',cindex2,fmt='%d')
np.savetxt('cindex3.txt',cindex3,fmt='%d');np.savetxt('cindex4.txt',cindex4,fmt='%d')
np.savetxt('cindex5.txt',cindex5,fmt='%d');np.savetxt('cindex6.txt',cindex6,fmt='%d')

#output overlaparea
np.savetxt('overlaparea.txt',overlaparea,fmt='%d')

#output b2overlap  
np.savetxt('b2overlap1.txt',b2overlap1,fmt='%d');np.savetxt('b2overlap2.txt',b2overlap2,fmt='%d')
np.savetxt('b2overlap3.txt',b2overlap3,fmt='%d');np.savetxt('b2overlap4.txt',b2overlap4,fmt='%d')
np.savetxt('b2overlap5.txt',b2overlap5,fmt='%d');np.savetxt('b2overlap6.txt',b2overlap6,fmt='%d')

#output cidim and c2innerL
cc=np.array([cidim1,c2innerL1]).T;np.savetxt('cidim+c2innerLP1.txt',cc,fmt='%d')
cc=np.array([cidim2,c2innerL2]).T;np.savetxt('cidim+c2innerLP2.txt',cc,fmt='%d')
cc=np.array([cidim3,c2innerL3]).T;np.savetxt('cidim+c2innerLP3.txt',cc,fmt='%d')
cc=np.array([cidim4,c2innerL4]).T;np.savetxt('cidim+c2innerLP4.txt',cc,fmt='%d')
cc=np.array([cidim5,c2innerL5]).T;np.savetxt('cidim+c2innerLP5.txt',cc,fmt='%d')
cc=np.array([cidim6,c2innerL6]).T;np.savetxt('cidim+c2innerLP6.txt',cc,fmt='%d')

#output over2in R abd L
over2in=np.array([over2inL1,over2inR1]).T;np.savetxt('over2inP1.txt',over2in,fmt='%d')
over2in=np.array([over2inL2,over2inR2]).T;np.savetxt('over2inP2.txt',over2in,fmt='%d')
over2in=np.array([over2inL3,over2inR3]).T;np.savetxt('over2inP3.txt',over2in,fmt='%d')
over2in=np.array([over2inL4,over2inR4]).T;np.savetxt('over2inP4.txt',over2in,fmt='%d')
over2in=np.array([over2inL5,over2inR5]).T;np.savetxt('over2inP5.txt',over2in,fmt='%d')
over2in=np.array([over2inL6,over2inR6]).T;np.savetxt('over2inP6.txt',over2in,fmt='%d')

#output c2bR and in2bR
tobR=np.array([c2bR1,in2bR1]).T;np.savetxt('2bR1.txt',tobR,fmt='%d')
tobR=np.array([c2bR2,in2bR2]).T;np.savetxt('2bR2.txt',tobR,fmt='%d')
tobR=np.array([c2bR3,in2bR3]).T;np.savetxt('2bR3.txt',tobR,fmt='%d')
tobR=np.array([c2bR4,in2bR4]).T;np.savetxt('2bR4.txt',tobR,fmt='%d')
tobR=np.array([c2bR5,in2bR5]).T;np.savetxt('2bR5.txt',tobR,fmt='%d')
tobR=np.array([c2bR6,in2bR6]).T;np.savetxt('2bR6.txt',tobR,fmt='%d')

