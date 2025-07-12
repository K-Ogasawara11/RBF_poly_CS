import numpy as np
from numpy import sin,cos,pi
from numba import jit
from scipy.spatial import KDTree
from numpy import pi,sin,cos,tan
import scipy
from params import ea



def knn_md(md):
    import numpy as np
    lad=np.loadtxt('../nodes/knn_md_N{}.txt'.format(md))
    lon=lad[:,0];lat=lad[:,1];x=lad[:,2];y=lad[:,3];z=lad[:,4];wgt=lad[:,5]
    return lon,lat,x,y,z,wgt

def knn_icos_nomod(md):
    import numpy as np
    lad=np.loadtxt('../nodes/knn_icos_nomodify_N{}.txt'.format(md))
    lon=lad[:,0];lat=lad[:,1];x=lad[:,2];y=lad[:,3];z=lad[:,4]
    return lon,lat,x,y,z


def bauer(n):
    import numpy as np
    import math
    m = np.floor(np.sqrt(n * math.pi) / 2) * 2 
    colat = np.arccos(1 - (2 * np.arange(0, n) + 1) / n)
    lon = np.mod(m * colat, math.tau)
    lat = math.tau / 4 - colat
    lon -=np.pi
    x=np.cos(lat)*np.cos(lon)
    y=np.cos(lat)*np.sin(lon)
    z=np.sin(lat)
    return lon, lat, x, y, z


def partition_xyz(x,y,z,lat,lon,N):
    import numpy as np
    from numpy import sin,cos,pi,tan
    pi4=pi/4
    pi34=pi/4*3
    a=1/np.sqrt(3)

    dim1=[]
    idim1=dim1.copy()
    idim2=dim1.copy()
    idim3=dim1.copy()
    idim4=dim1.copy()
    idim5=dim1.copy()
    idim6=dim1.copy()
    for i in range(N):
        if lon[i]<=pi4 and lon[i]>-pi4:
            Y=a*(z[i]/x[i])
            if x[i]==0:
                if lat[i]>pi4:
                    idim5=np.append(idim5,i)
                elif lat[i]<=-pi4:
                    idim6=np.append(idim6,i)
            elif Y>=a:
                idim5=np.append(idim5,i)
            elif Y<=-a:
                idim6=np.append(idim6,i)
            else:
                idim1=np.append(idim1,i)
        elif lon[i]<=pi34 and lon[i]>pi4:
            Y=a*(z[i]/y[i])
            if Y>=a:
                idim5=np.append(idim5,i)
            elif Y<=-a:
                idim6=np.append(idim6,i)
            elif z[i]==1:
                idim5=np.append(idim5,i)
            elif z[i]==-1:
                idim6=np.append(idim6,i)
            else:
                idim2=np.append(idim2,i)
        elif lon[i]<=-pi34:
            
            Y=a*(-z[i]/x[i])
            if Y>=a:
                idim5=np.append(idim5,i)
                
            elif Y<=-a:
                idim6=np.append(idim6,i)
            elif z[i]==1:
                idim5=np.append(idim5,i)
            elif z[i]==-1:
                idim6=np.append(idim6,i)
            elif Y<a and Y>-a:
                idim3=np.append(idim3,i)
                
        elif lon[i]>pi34:
            Y=a*(-z[i]/x[i])
            if Y>=a:
                idim5=np.append(idim5,i)
            elif Y<=-a:
                idim6=np.append(idim6,i)
            elif z[i]==1:
                idim5=np.append(idim5,i)
            elif z[i]==-1:
                idim6=np.append(idim6,i)
            elif Y<a and Y>-a:
                idim3=np.append(idim3,i)
                
        elif lon[i]<=-pi4 and lon[i]>-pi34:
            
            Y=a*(-z[i]/y[i])
            if Y>=a:
                idim5=np.append(idim5,i)
            elif Y<=-a:
                idim6=np.append(idim6,i)
            elif z[i]==1:
                idim5=np.append(idim5,i)
            elif z[i]==-1:
                idim6=np.append(idim6,i)
            else:
                idim4=np.append(idim4,i)
    
    dim=np.array([])
    dim=np.append(dim,idim1)
    dim=np.append(dim,idim2)
    dim=np.append(dim,idim3)
    dim=np.append(dim,idim4)
    dim=np.append(dim,idim5)
    dim=np.append(dim,idim6)
    print('node＝＝sum(idim.size)',dim.size==N,dim.size,N)
    
    return np.int64(idim1),np.int64(idim2),np.int64(idim3),np.int64(idim4),np.int64(idim5),np.int64(idim6)




def ca2cs_eq(u,v,lat,lon,P):
    
    a=ea/np.sqrt(3)
    if P==2:
        clon=lon.copy()
        clon=lon-pi/2
    elif P==3:
        clon=lon.copy()
        #clon=lon-pi
        #print('P3 max clon',np.amin(clon/pi*180))
        for i in range(lon.shape[0]):        
            if lon[i]<0:
                clon[i]=clon[i]+pi*2
        clon-=pi    
    elif P==4:
        clon=lon.copy()
        clon=lon+pi/2
        #print('P4 max clon',np.amax(clon/pi*180))
    else:
        clon=lon
   
    u1=a*(u/(cos(lat)*cos(clon)**2))
    v1=a*(u*tan(clon)*tan(lat)+v/cos(lat))/(cos(lat)*cos(clon))
    
    return u1,v1

def ca2cs_eqv2(u,v,coslat,coslon,tanlat,tanlon):
    
    a=ea/np.sqrt(3)
    u1=a*(u/(coslat*coslon**2))
    v1=a*(u*tanlon*tanlat+v/coslat)/(coslat*coslon)
    return u1,v1


def calc_covariant_eqv2(u1,v1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    cu1=1/a**2*((clat**2*clon**2*(clon**2+slat**2*slon**2))*u1-clat**3*clon**2*slat*slon*v1)
    #cu1=-clat**3*clon**2*slat*slon*v1
    cv1=1/a**2*((-clat**3*clon**2*slat*slon*u1)+clon**2*clat**4*v1)
    return cu1,cv1 


def co2contra_eqv2(cu1,cv1,clat,clon,slat,slon):
    
    a=ea/np.sqrt(3)

    #clat=cos(lat);clon=cos(cslon);slon=sin(cslon);slat=sin(lat)    
    u1=a**2*(clon**2*clat**4*cu1+clat**3*clon**2*slat*slon*cv1)/(clat**6*clon**6)
    v1=a**2*((clat**3*clon**2*slat*slon*cu1)+(clat**2*clon**2*(clon**2+slat**2*slon**2))*cv1)/(clat**6*clon**6)
    #print(u1.shape,v1.shape)
    return u1,v1


def cs2ca_eqv2(u1,v1,coslat,coslon,sinlat,sinlon):
    
    a=ea/np.sqrt(3)
    u=(u1*coslat*coslon**2)/a
    v=(-u1*sinlat*sinlon*coslat*coslon+v1*coslon*coslat**2)/a
    return u,v

def ca2cs_npv2(u,v,coslat,coslon,sinlat,sinlon):
    
    a=ea/np.sqrt(3)
    u1=(u*sinlat*coslon-v*sinlon)/sinlat**2*a
    v1=(u*sinlat*sinlon+v*coslon)/sinlat**2*a
    return u1,v1

def calc_covariant_npv2(u1,v1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)    
    cu1=(slat**2*(clon**2+slat**2*slon**2)*u1+clat**2*slat**2*slon*clon*v1)/a**2
    cv1=(clat**2*slat**2*slon*clon*u1+slat**2*(slon**2+clon**2*slat**2)*v1)/a**2
    return cu1,cv1

def co2contra_npv2(cu1,cv1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    u1=a**2*(slat**2*(slon**2+clon**2*slat**2)*cu1-clat**2*slat**2*slon*clon*cv1)/(slat**6)
    v1=a**2*(-clat**2*slat**2*slon*clon*cu1+slat**2*(clon**2+slat**2*slon**2)*cv1)/(slat**6)
    return u1,v1

def ca2cs_np_eang(u,v,lat,lon,x,y,z):
   
    alpha=np.arctan((y/z))
    beta=np.arctan((-x/z))
    u1=cos(alpha)**2*(u*sin(lat)*cos(lon)-v*sin(lon))/sin(lat)**2
    v1=cos(beta)**2*(u*sin(lat)*sin(lon)+v*cos(lon))/sin(lat)**2
    return u1,v1

def cs2ca_npv2(u1,v1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    u=slat*(u1*clon+v1*slon)/a
    v=slat*(-u1*slat*slon+v1*slat*clon)/a
    return u,v

def calc_covariant_spv2(u1,v1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    cu1=(slat**2*(clon**2+slat**2*slon**2)*u1-clat**2*slat**2*slon*clon*v1)/a**2
    cv1=(-clat**2*slat**2*slon*clon*u1+slat**2*(slon**2+clon**2*slat**2)*v1)/a**2
    return cu1,cv1

def co2contra_spv2(cu1,cv1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    u1=a**2*(slat**2*(slon**2+clon**2*slat**2)*cu1+clat**2*slat**2*slon*clon*cv1)/(slat**6)
    v1=a**2*(clat**2*slat**2*slon*clon*cu1+slat**2*(clon**2+slat**2*slon**2)*cv1)/(slat**6)
    return u1,v1

def ca2cs_spv2(u,v,clat,clon,slat,slon):
    
    a=ea/np.sqrt(3)
    u1=(-u*slat*clon+v*slon)/slat**2*a
    v1=(u*slat*slon+v*clon)/slat**2*a
    return u1,v1


def ca2cs_sp_eang(u,v,lat,lon,x,y,z):
    
    alpha=np.arctan((-y/z))
    beta=np.arctan((-x/z))
    u1=cos(alpha)**2*(-u*sin(lat)*cos(lon)+v*sin(lon))/sin(lat)**2
    v1=cos(beta)**2*(u*sin(lat)*sin(lon)+v*cos(lon))/sin(lat)**2
    return u1,v1

def cs2ca_spv2(u1,v1,clat,clon,slat,slon):
    a=ea/np.sqrt(3)
    u=slat/a*(-u1*clon+v1*slon)
    v=slat/a*(u1*slat*slon+v1*slat*clon)
    return u,v


def find_corner(x,y):
    r=x**2+y**2
    ind=(np.where(r==np.amax(r)))[0]
    return ind[0]

def plane_ordering(index,x,y):
    corner=find_corner(x,y)
    

if __name__ == '__main__':
    print('work')