from math import tau, pi
import numpy as np
from numpy import cos, sin, tan, arcsin, arccos, exp
import sys
from params import ea,day_in_sec,grav,omega,alpha,rea


def cosbell2(lon, lat, lon0, lat0, su0, alpha , t, a, R):
    
    if t==0:
        lonc = lon0
        latc = lat0
        print('time=0')
    else:
        lonc, latc = rotate2(lon0+ su0*t , lat0, alpha)
    #print('経度緯度は',lonc,latc)
    h0 = 1
    r = arccos(sin(latc) * sin(lat) + cos(latc) * cos(lat) * cos(lon - lonc))
    #print(r.shape)
    h = h0 / 2 * (1 + cos(pi * r / R))
    h[r >= R] = 0
    
    return h

def gaussbell2(lon, lat, lon0, lat0, su0, alpha , t, a, R):
    
    if t==0:
        lonc = lon0
        latc = lat0
        print('time=0')
    else:
        lonc, latc = rotate2(lon0+ su0*t , lat0, alpha)
    #print('経度緯度は',lonc,latc)
    h0 = 1000.0
    r = arccos(sin(latc) * sin(lat) + cos(latc) * cos(lat) * cos(lon - lonc))
    #print(r.shape)
    h = h0 *np.exp(-(2.25*r/R)**2)
    #h = h0 *np.exp(-6*r**2)
    
    return h


def gaussbell3(x,lon0,lon,lat):
    lonc=lon-lon0
    for i in range(lon.size):
       if lonc[i]<-pi:
            lonc[i]+=pi*2
    x2=cos(lonc)*cos(lat)
    rho =np.arccos(x2)
    h=np.exp(-6*rho**2)
    return h
def rotate2(lon, lat, alpha):
    if alpha==0:
        rolon=lon
        rolat=lat
    else:
        x=cos(lat)*cos(lon)
        y=cos(lat)*sin(lon)
        z=sin(lat)
        xb=-sin(alpha)*z+cos(alpha)*x
        zb=cos(alpha)*z+sin(alhpa)*x
        yb=y
        rolon=np.arctan2(yb/xb)
        rolat=arcsin(zb)

    return rolon, rolat


def rotate(lon, lat, alpha):
    if alpha == 0:
        rotlon = lon
        rotlat = lat
    else:
        test = sin(lat) * cos(alpha) - cos(lat) * cos(lon) * sin(alpha)
        if test > 1:
            rotlat = pi / 2
        elif test < -1:
            rotlat = -pi / 2
        else:
            rotlat = arcsin(test)
        test = cos(rotlat)
        if test == 0:
            rotlon = 0.0
        else:
            test = sin(lon) * cos(lat) / test
            if test > 1:
                rotlon = pi / 2
            elif test < -1:
                rotlon = -pi / 2
            else:
                rotlon = arcsin(test)
        test = cos(alpha) * cos(lon) * cos(lat) + sin(alpha) * sin(lat)
        if test < 0:
            rotlon = pi - rotlon
        rotlon = (rotlon + tau) % tau
    return rotlon, rotlat


    

def cosbell(lon, lat, lon0=0.0, lat0=0.0, su0=0.0, alpha=0 , t=0.0, a=1.0, R=1/3):
    if t != 0.0:
        lonc, latc = rotate(lon0 + su0 * cos(alpha) * t / a, lat0, alpha)
    else:
        lonc = lon0
        latc = lat0
    #lonc -= tau / 4
    #print('経度緯度は',lonc,latc)
    h0 = 1000.0
    r = arccos(sin(latc) * sin(lat) + cos(latc) * cos(lat) * cos(lon - lonc))
    h = h0 / 2 * (1 + cos(pi * r / R))
    h[r >= R] = 0
    return h

def case2(lon,lat):
    iea=6.37122e6
    u0=2*pi*iea/(12*day_in_sec)
    h0=2.94e4
    u=u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
    v=-u0*sin(lon)*sin(alpha)
    h=h0-(iea*omega*u0+u0**2/2)*(-cos(lon)*cos(lat)*sin(alpha)+sin(lat)*cos(alpha))**2

    return u,v,h/grav

def case5(lon,lat,mount):
    iea=6.37122e6
    u0=20
    if mount=='piece':
        h0=5960 
    elif mount=='gauss':
        h0=5400
    u=u0*(np.cos(lat))
    v=np.zeros(lat.size)#-u0*sin(lon)*sin(alpha)
    h=h0-(iea*omega*u0+u0**2/2)*(-cos(lon)*cos(lat)*sin(alpha)+sin(lat)*cos(alpha))**2/grav

    return u,v,h


def case5_surf(lon,lat,mount):
    hs0=2000
    latc=pi/6
    lonc=3*pi/2
    rr=pi/9
    lon2=lon.copy()
    for i in range(lon.size):
        if lon[i]<0:
            lon2[i]+=pi*2
    r2=np.zeros(lat.size)
    for i in range(lat.size):
        r2[i]=min(rr**2,(lon2[i]-lonc)**2+(lat[i]-latc)**2)
    r=np.sqrt(r2)
    if mount=='piece':
        hs=hs0*(1-r/rr)
    elif mount=='gauss':
        hs=hs0*np.exp(-2.5*(r/rr)**2)
    return hs

def case6(lon,lat):
    omega6=7.848e-6
    K=omega6
    R=4
    h0=8e3
    u=rea*omega6*cos(lat)+rea*K*cos(lat)**(R-1)*(R*sin(lat)**2-cos(lat)**2)*cos(R*lon)
    v=-rea*K*R*cos(lat)**(R-1)*sin(lat)*sin(R*lon)

    Agh=omega6*0.5*(2*omega+omega6)*cos(lat)**2\
        +0.25*K**2*cos(lat)**(2*R)\
            *((R+1)*cos(lat)**2+((2*R**2)-R-2)-2*R**2/cos(lat)**(2))
    Bgh=2*(omega+omega6)*K/((R+1)*(R+2))\
        *cos(lat)**(R)\
            *((R**2+2*R+2)-(R+1)**2*cos(lat)**2)
    Cgh=0.25*K**2*cos(lat)**(2*R)*\
        ((R+1)*cos(lat)**2-(R+2))

    gh=grav*h0+rea**2*Agh+rea**2*Bgh*cos(R*lon)+rea**2*Cgh*cos(2*R*lon)

    return u,v,gh/grav

def galewsky(lon,lat):
    fn=lat.size
    umax=80
    lat0=pi/7
    lat1=pi/2-lat0
    en=np.exp(-4/(lat1-lat0)**2)
    v=np.zeros(fn)
    u=np.zeros(fn)
    for i in range(fn):
        if lat[i]<=lat0:
            u[i]=0
        elif lat[i] >=lat1:
            u[i]=0
        else: 
            u[i]=umax/en*np.exp(1/((lat[i]-lat0)*(lat[i]-lat1)))
    h=np.loadtxt('galewsky_init_hightN{}.txt'.format(fn))
    #data=np.loadtxt('init_galewsky_helpixN{}.txt'.format(fn))
    #h=data[:,0]
    #u=data[:,1]

    return u,v,h

def galewskyb(N):
    deta=np.loadtxt('init_galewsky_baranse_helpixN{}.txt'.format(N))
    h=deta[0]
    u=deta[1]
    v=np.zeros(N)
    return h,u,v

def calc_h(lat):
    omega=7.292e-5;g=9.80616;h0=10000
    umax=80
    a=6.37122e6
    lat0=pi/7;lat1=pi/2-lat0;en=np.exp(-4/(lat1-lat0)**2)
    u=np.zeros(lat.size)
    for i in range(lat.size):
        if lat[i]>lat0 and lat[i]<lat1: 
            u[i]=umax/en*np.exp(1/((lat[i]-lat0)*(lat[i]-lat1)))
        if lat[i]<=lat0:
            u[i]=0
        elif lat[i]>lat1:
            u[i]=0
            
    ghd=a*u*2*omega*sin(lat)+tan(lat)*u**2
    return ghd

if __name__ == '__main__':
    import numpy as np
    print('work')