import numpy as np
import md_cube as cube
from scipy.spatial import KDTree

def calc_x1y1(xb,yb,zb,pn):
    a=1/np.sqrt(3)
    if pn==1:
        X1=a*(yb/xb)
        Y1=a*(zb/xb)
    elif pn==2:
        X1=a*(-xb/yb)
        Y1=a*(zb/yb)
    elif pn==3:
        X1=a*(yb/xb)
        Y1=a*(-zb/xb)
    elif pn==4:
        X1=a*(-xb/yb)
        Y1=a*(-zb/yb)
    elif pn==5:
        Y1=a*(-xb/zb)
        X1=a*(yb/zb)
    elif pn==6:
        Y1=a*(-xb/zb)
        X1=a*(-yb/zb)
    return X1,Y1

def convert_stensil_index(x,y,z,x1,y1,z1,idim1,pn,n):
    xc1,yc1=calc_x1y1(x1,y1,z1,pn)
    corner=cube.find_corner(xc1,yc1)
    xc=xc1[corner];yc=yc1[corner]
    index=np.zeros((idim1.size,n))
    tree = KDTree(np.array([x,y,z]).T)
    for i in range(idim1.size):
        value,indexb=tree.query([x1[i],y1[i],z1[i]],k=n)
        index[i,:]=np.array(indexb)
    uni_index=np.int64(np.unique(index))
    xb=x[uni_index];yb=y[uni_index];zb=z[uni_index]
    xcb1,ycb1=calc_x1y1(xb,yb,zb,pn)
    tree = KDTree(np.array([xcb1,ycb1]).T)
    value,indexb=tree.query([xc,yc],k=uni_index.size)
    cindexb=np.int64(np.array(indexb))
    cindex=uni_index[cindexb]
    return cindex

def calc_overlaparea(cindex1,cindex2,cindex3,cindex4,cindex5,cindex6):
    overlaparea=[]
    both=np.intersect1d(cindex1,cindex2);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex1,cindex4);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex1,cindex5);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex1,cindex6);overlaparea=np.append(overlaparea,both)

    both=np.intersect1d(cindex2,cindex3);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex2,cindex5);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex2,cindex6);overlaparea=np.append(overlaparea,both)

    both=np.intersect1d(cindex3,cindex4);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex3,cindex5);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex3,cindex6);overlaparea=np.append(overlaparea,both)

    both=np.intersect1d(cindex4,cindex5);overlaparea=np.append(overlaparea,both)
    both=np.intersect1d(cindex4,cindex6);overlaparea=np.append(overlaparea,both)
    
    overlapareab=np.unique(overlaparea)
    return np.int64(overlapareab)

def ordering_idim(x1,y1,z1,idim1,pn):
    xc1,yc1=calc_x1y1(x1,y1,z1,pn)
    corner=cube.find_corner(xc1,yc1)
    xc=xc1[corner];yc=yc1[corner]
    tree = KDTree(np.array([xc1,yc1]).T)
    value,indexb=tree.query([xc,yc],k=idim1.size)
    cidimb=np.int64(np.array(indexb))
    cidim=idim1[cidimb]
    return np.int64(cidim)

def calc_initeraction_index(overlaparea,cidim1,cindex1):
    interaction_area1=np.intersect1d(cidim1,overlaparea)
    over2inb1=np.intersect1d(cindex1,overlaparea)
    over2inL1=over2inb1.copy()
    over2inR1=over2inb1.copy()
    for i in range(over2inb1.size):
        over2inL1[i]=int(np.where(cindex1==over2inb1[i])[0])
        over2inR1[i]=int(np.where(overlaparea==over2inb1[i])[0])
    c2innerb1=np.zeros(cidim1.size)
    for i in range(cidim1.size):
        c2innerb1[i]=np.where(cindex1==cidim1[i])[0]
    c2innerL1=np.int64(c2innerb1)

    cs2b_indexb1=np.zeros(interaction_area1.size)
    for i in range(interaction_area1.size):
        cs2b_indexb1[i]=np.where(cidim1==interaction_area1[i])[0]
    c2bR=np.int64(cs2b_indexb1)

    for i in range(interaction_area1.size):
        cs2b_indexb1[i]=np.where(overlaparea==interaction_area1[i])[0]
    b2overlap=np.int64(cs2b_indexb1)

    in2bR=np.zeros(c2bR.size)
    for i in range(c2bR.size):
        
    	in2bR[i]=int(c2innerL1[c2bR[i]])
    in2bR=np.int64(in2bR)
    
    in2cR=np.zeros(cidim1.size)
    for I in range(cidim1.size):
    	in2cR[i]=np.where(cindex1==cidim1[i])[0]
    in2cR=np.int64(in2cR)

    return np.int64(over2inL1),np.int64(over2inR1),np.int64(c2innerL1),np.int64(c2bR),np.int64(b2overlap),in2bR,in2cR

def calc_overlaparea2(cindex1,cidim1,cindex2,cidim2,cindex3,cidim3,cindex4,cidim4,cindex5,cidim5,cindex6,cidim6):

    overlaparea=[]
    uni=np.unique(cindex1)
    for i in range(cidim1.size):
        ind=(np.where(uni==cidim1[i]))[0]
        uni[ind]=-1
    uni1=np.unique(uni).copy()
    
    uni=np.unique(cindex2)
    for i in range(cidim2.size):
        ind=(np.where(uni==cidim2[i]))[0]
        uni[ind]=-1
    uni2=np.unique(uni).copy()

    uni=np.unique(cindex3)
    for i in range(cidim3.size):
        ind=(np.where(uni==cidim3[i]))[0]
        uni[ind]=-1
    uni3=np.unique(uni).copy()

    uni=np.unique(cindex4)
    for i in range(cidim4.size):
        ind=(np.where(uni==cidim4[i]))[0]
        uni[ind]=-1
    uni4=np.unique(uni).copy()

    uni=np.unique(cindex5)
    for i in range(cidim5.size):
        ind=(np.where(uni==cidim5[i]))[0]
        uni[ind]=-1
    uni5=np.unique(uni).copy()

    uni=np.unique(cindex6)
    for i in range(cidim6.size):
        ind=(np.where(uni==cidim6[i]))[0]
        uni[ind]=-1
    uni6=np.unique(uni).copy()
    
    return uni1,uni2,uni3,uni4,uni5,uni6

def calc_initeraction_indexv2(overlaparea,cidim1,cindex1,uni1):
    interaction_area1=np.intersect1d(uni1,overlaparea)
    over2inL1=interaction_area1.copy()
    over2inR1=interaction_area1.copy()
    for i in range(interaction_area1.size):
        over2inL1[i]=int(np.where(cindex1==interaction_area1[i])[0])
        over2inR1[i]=int(np.where(overlaparea==interaction_area1[i])[0])
    return np.int64(over2inL1),np.int64(over2inR1)



if __name__ == '__main__':
    print('work')