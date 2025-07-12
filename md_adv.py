import numpy as np
from numba import jit
from scipy.spatial import KDTree
from numpy import pi,sin,cos,tan
import scipy



def calc_R_phs(sx,x,y,z,index,n,pn):
    a=1/np.sqrt(3)
    R=np.zeros(((sx.size,n,n)))
    dx=np.zeros((sx.size,n))
    dy=np.zeros((sx.size,n))
    dz=np.zeros((sx.size,n))
    stx=np.zeros((sx.size,n))
    sty=np.zeros((sx.size,n))
    stx0=np.zeros((sx.size,n))
    sty0=np.zeros((sx.size,n))
    if pn==1:
        px=a*(y/x);py=a*(z/x)
    elif pn==2:
        px=a*(-x/y);py=a*(z/y)
    elif pn==3:
        px=a*(y/x);py=a*(-z/x)
    elif pn==4:
        px=a*(-x/y);py=a*(-z/y)
    elif pn==5:
        px=a*(y/z);py=a*(-x/z)
    elif pn==6:
        px=a*(-y/z);py=a*(-x/z)
    for i in range(sx.size):
        for j in range(n):
            inda=int(index[i,j])
            for k in range(n):
                indb=int(index[i,k])
                R[i,j,k]=(x[inda]-x[indb])**2+(y[inda]-y[indb])**2+(z[inda]-z[indb])**2
    for i in range(sx.size):
        for j in range(n):
            indb=int(index[i,j])
            inda=int(index[i,0])
            dx[i,j]=x[inda]-x[indb]
            dy[i,j]=y[inda]-y[indb]
            dz[i,j]=z[inda]-z[indb]
            stx[i,j]=px[inda]
            sty[i,j]=py[inda]
            stx0[i,j]=px[indb]
            sty0[i,j]=py[indb]
    r=np.sqrt(a**2+stx**2+sty**2)
    return R,dx,dy,dz,stx,sty,stx0,sty0,r

def calc_poly(px,py,deg,n,projection):
    a=1/np.sqrt(3)
    pol=int((deg+1)*(deg+2)/2)
    poly=np.ones((n,pol))
    alpha=np.arctan(px[0]/a)
    beta=np.arctan(py[0]/a)
    dxpoly=np.zeros(pol)
    dypoly=np.zeros(pol)
    a=1/np.sqrt(3)
    stx=px-px[0]
    sty=py-py[0]
    X=np.vstack((stx,stx))
    Y=np.vstack((sty,sty))
    for i in range(pol):
        X=np.vstack((X,stx))
        Y=np.vstack((Y,sty))
    X[0,:]=1
    Y[0,:]=1
    Y=np.cumprod(Y,axis=0)
    X=np.cumprod(X,axis=0)
    Y=Y.T;X=X.T
    col=0
    py=np.zeros(n)
    for j in range(deg+1):
        for k in range(j+1):
            py=Y[:,j-k]
            poly[:,col]=py*X[:,k]
            col+=1
    if deg>1:
        if projection=='edist':
            dxpoly[2]=1
            dypoly[1]=1
        elif projection=='eang':
            dxpoly[2]=a/cos(alpha)**2
            dypoly[1]=a/cos(beta)**2
    return poly,dxpoly,dypoly


def calc_ARA_phs(R,N,n,stx,sty,deg,projection):
    pol=int((deg+1)*(deg+2)/2)
    R1=np.sqrt(R)
    A=R1**3
    RA=np.zeros(((N,n+pol,n+pol)))
    dxpols=np.zeros((N,n+pol))
    dypols=np.zeros((N,n+pol))
    RA[:,0:n,0:n]=A
    for i in range(N):
        poly,dxpoly,dypoly=calc_poly(stx[i,:],sty[i,:],deg,n,projection)
        RA[i,n:n+pol,0:n]=poly.T
        RA[i,0:n,n:n+pol]=poly
        dxpols[i,n:n+pol]=dxpoly
        dypols[i,n:n+pol]=dypoly
    return A,RA,dxpols,dypols

def calc_coefficient(a,x,y,r,pn):
    r3=r**3;r2=r**2
    x2=x**2;y2=y**2
    if pn==1:
        xx=-a*x/r3;xy=(r2-x2)/r3;xz=-y*x/r3
        yx=-a*y/r3;yy=-y*x/r3;yz=(r2-y2)/r3   
    elif pn==2:
        xx=(-r2+x2)/r3;xy=-a*x/r3;xz=-y*x/r3
        yx=x*y/r3;yy=-a*y/r3;yz=(r2-y2)/r3 
    elif pn==3:
        xx=a*x/r3;xy=(-r2+x2)/r3;xz=-x*y/r3
        yx=a*y/r3;yy=(y*x)/r3;yz=(r2-y2)/r3
    elif pn==4:
        xx=(r2-x2)/r3;xy=a*x/r3;xz=-y*x/r3
        yx=-x*y/r3;yy=a*y/r3;yz=(r2-y2)/r3
    elif pn==5:
        xx=x*y/r3;xy=(r2-x2)/r3;xz=-a*x/r3
        yx=(-r2+y2)/r3;yy=-x*y/r3;yz=-a*y/r3
    elif pn==6:
        xx=(-x*y)/r3;xy=(r2-x2)/r3;xz=a*x/r3
        yx=(r2-y2)/r3;yy=-x*y/r3;yz=a*y/r3
    return xx,xy,xz,yx,yy,yz


def calc_Dmat_poly(RA,R,r,dx,dy,dz,stx,sty,bdx,bdy,sx,n,P):
    a=1/np.sqrt(3)
    Dx=np.zeros((sx.size,n))
    Dy=np.zeros((sx.size,n))
    xx,xy,xz,yx,yy,yz=calc_coefficient(a,stx,sty,r,P)
    for i in range(sx.size):
        dphir=3*R[i,0,:]**(1/2) 
        bdxb=dphir*(dx[i,:]*xx[i,:]+dy[i,:]*xy[i,:]+dz[i,:]*xz[i,:])
        bdyb=dphir*(dx[i,:]*yx[i,:]+dy[i,:]*yy[i,:]+dz[i,:]*yz[i,:])
        bdx[i,0:n]=bdxb
        bdy[i,0:n]=bdyb
        
        dxb=np.linalg.solve(RA[i,:,:],bdx[i,:])
        dyb=np.linalg.solve(RA[i,:,:],bdy[i,:])
        Dx[i,:]=dxb[0:n];Dy[i,:]=dyb[0:n]
    return Dx,Dy



def calc_Dmatrix_phs(x,y,z,x1,y1,z1,idim,pn,n,deg,projection):
    index=np.zeros((idim.size,n))
    tree = KDTree(np.array([x,y,z]).T)
    for i in range(idim.size):
        value,indexb=tree.query([x1[i],y1[i],z1[i]],k=n)
        index[i,:]=np.array(indexb)
    R,dx,dy,dz,stx,sty,stx0,sty0,r=calc_R_phs(idim,x,y,z,index,n,pn)
    A,RA,dxpols,dypols=calc_ARA_phs(R,idim.size,n,stx0,sty0,deg,projection) 
    if projection=='edist':
        Dx,Dy=calc_Dmat_poly(RA,R,r,dx,dy,dz,stx,sty,dxpols,dypols,idim,n,pn)
    elif projection=='eang':
        Dx,Dy=calc_Dmat_poly_eang(RA,R,r,dx,dy,dz,stx,sty,dxpols,dypols,idim,n,pn)
    

    return Dx,Dy,r,index


def calc_dh_adv(SDx,SDy,u,v,p_vec2all,h):
    pdh=(SDx@h)*u
    pdh+=(SDy@h)*v
    dh=p_vec2all@pdh
    return dh

def make_transform_matrix(idim,N):
    sn=idim.size
    vals=np.ones(sn)
    row=np.zeros(sn)
    col=np.zeros(sn)
    
    for i in range(sn):
        row[i]=i
        col[i]=idim[i]
    p_vec2all=scipy.sparse.csr_matrix((vals,(col,row)),shape = (N,sn))
    row=np.zeros(sn)
    col=np.zeros(sn)
    for i in range(sn):
        row[i]=idim[i]
        col[i]=i
    p_all2vec=scipy.sparse.csr_matrix((vals,(col,row)),shape = (sn,N))

    return p_vec2all,p_all2vec

@jit(nopython=True)
def make_row_col(index,val,N,n):
    row=np.zeros(N*n)
    col=np.zeros(N*n)
    vals=np.zeros(N*n)
    count=0
    for i in range(N):
        for j in range(n):
            ind=int(index[i,j])
            row[count]=i
            col[count]=ind
            vals[count]=val[i,j]
            count+=1
    return row,col,vals

def convert_sparse(index,Dx,Dy,idim,N,n):
    row,col,vals=make_row_col(index,Dx,idim.size,n)
    SDx=scipy.sparse.csr_matrix((vals,(row,col)),shape = (idim.size,N))
    row,col,vals=make_row_col(index,Dy,idim.size,n)
    SDy=scipy.sparse.csr_matrix((vals,(row,col)),shape = (idim.size,N))
    return SDx,SDy

@jit(nopython=True)
def make_row_col_pc(index,vals,n,N):
    val2=vals.reshape(N*n)
    col=index.reshape(N*n)
    count=0
    row=np.zeros(N*n)
    for i in range(N):
        for j in range(n):
            row[count]=i
            count+=1
    return row,col,val2

def calc_index_pc(N,n,x,y,z):
    index=np.zeros((N,n))
    tree = KDTree(np.array([x,y,z]).T)
    for i in range(N):
        value,indexb=tree.query([x[i],y[i],z[i]],k=n)
        index[i,:]=np.array(indexb)
    return index    

@jit(nopython=True)
def calc_R_pc(N,n,x,y,z,indexb):
    R=np.zeros(((N,n,n)))
    for i in range(N):
        for j in range(n):
            for l in range(n):
                a=int(indexb[i,j]);b=int(indexb[i,l])
                R[i,j,l]=(x[a]-x[b])**2+(y[a]-y[b])**2+(z[a]-z[b])**2
    return R

# Fornberg and Lehto 2011
def calc_hiv_FH2011(d,k,ep,R,N,n):
    import math
    import numpy as np
    import scipy.special as sp
    laplas2=np.zeros((N,n))
    ep2=ep**2
    LA=np.zeros((n+1,n+1))
    BB=np.zeros(n+1)
    LA[n,:]=1
    LA[:,n]=1
    LA[n,n]=0
    for i in range(N):
        ra=R[i,:,:]
        rb=R[i,0,:]
        p=calc_p(rb,ep2,k,n,d)
        laplas=ep2**k*p
        BB[0:n]=(laplas*np.exp(-ep2*rb)).T
        LA[0:n,0:n]=np.exp(-ep2*ra)
        lap=np.linalg.solve(LA,BB).T#;print(lap.shape)
        laplas2[i,:]=lap[0:n]
    return laplas2

# Fornberg and Lehto 2011
def calc_p(r,ep,k,n,d):
    sr=r.size
    pk=np.zeros((n+1,sr))
    epr2=(ep*r)
    pk[0]=1
    pk[1]=4*(epr2)-2*d
    for i in range(2,k+1):
        pk[i]=4*(epr2-2*(i-1)-d/2)*pk[i-1]-8*(i-1)*(2*(i-1)-2+d)*pk[i-2]
    p=pk[k].copy()
    return p


def ca2cs_eq(u,v,lat,lon,P):
    ea=1
    a=ea/np.sqrt(3)
    if P==2:
        clon=lon-pi/2
    elif P==3:
        clon=lon.copy()
        for i in range(lon.size):
            if lon[i]<0:
                clon[i]=lon[i]+pi*2
        clon-=pi    
    elif P==4:
        clon=lon+pi/2
        
    else:
        clon=lon
   
    u1=a*(u/(cos(lat)*cos(clon)**2))
    v1=a*(u*tan(clon)*tan(lat)+v/cos(lat))/(cos(lat)*cos(clon))
    return u1,v1


def ca2cs_np(u,v,lat,lon):
    ea=1
    a=ea/np.sqrt(3)
    u1=(u*sin(lat)*cos(lon)-v*sin(lon))/sin(lat)**2*a
    v1=(u*sin(lat)*sin(lon)+v*cos(lon))/sin(lat)**2*a
    return u1,v1

def ca2cs_sp(u,v,lat,lon):
    ea=1
    a=ea/np.sqrt(3)
    u1=(-u*sin(lat)*cos(lon)+v*sin(lon))/sin(lat)**2*a
    v1=(u*sin(lat)*sin(lon)+v*cos(lon))/sin(lat)**2*a
    return u1,v1

def ca2cs_eq_eig(u,v,lat,lon,P,n):
    ea=1
    a=ea/np.sqrt(3)
    if P==2:
        clon=lon-pi/2
    elif P==3:
        clon=lon.copy()
        for i in range(lon.shape[0]):
            for j in range(n):
                if lon[i,j]<0:
                    clon[i,j]=lon[i,j]+pi*2
        clon-=pi    
    elif P==4:
        clon=lon+pi/2
    else:
        clon=lon
   
    u1=a*(u/(cos(lat)*cos(clon)**2))
    v1=a*(u*tan(clon)*tan(lat)+v/cos(lat))/(cos(lat)*cos(clon))

    return u1,v1

def metric_eq_flux(Dx,Dy,index,idim1,u,v,lat,lon,n,P):
    ub=np.zeros((idim1.size,n))
    latb=np.zeros((idim1.size,n))
    lonb=np.zeros((idim1.size,n))
    vb=np.zeros((idim1.size,n))
    for i in range(idim1.size):
        for j in range(n):
            ind=int(index[i,j])
            ub[i,j]=u[ind]
            vb[i,j]=v[ind]
            latb[i,j]=lat[ind]
            lonb[i,j]=lon[ind]
    u1,v1=ca2cs_eq_eig(ub,vb,latb,lonb,P,n)
    a=1/np.sqrt(3)
    g1=calc_jacob(lonb,latb,P,n)
    Dx=Dx*g1*u1/a**2
    Dy=Dy*g1*v1/a**2
    return Dx,Dy


def metric_np_flux(Dx,Dy,index,idim1,u,v,lat,lon,n):
    P=5
    a=1/np.sqrt(3)
    ub=np.zeros((idim1.size,n))
    latb=np.zeros((idim1.size,n))
    lonb=np.zeros((idim1.size,n))
    vb=np.zeros((idim1.size,n))
    ga=np.ones((idim1.size,n))
    for i in range(idim1.size):
        for j in range(n):
            ind=int(index[i,j])
            ub[i,j]=u[ind]
            vb[i,j]=v[ind]
            latb[i,j]=lat[ind]
            lonb[i,j]=lon[ind]
    ga=calc_jacob(lonb,latb,P,n)
    u1,v1=ca2cs_np(ub,vb,latb,lonb)
  
    Dx=ga*u1*Dx/a**2
    Dy=ga*v1*Dy/a**2
    return Dx,Dy

def metric_sp_flux(Dx,Dy,index,idim1,u,v,lat,lon,n):
    P=6
    a=1/np.sqrt(3)
    ub=np.zeros((idim1.size,n))
    latb=np.zeros((idim1.size,n))
    lonb=np.zeros((idim1.size,n))
    vb=np.zeros((idim1.size,n))
    ga=np.ones((idim1.size,n))
    for i in range(idim1.size):
        for j in range(n):
            ind=int(index[i,j])
            ub[i,j]=u[ind]
            vb[i,j]=v[ind]
            latb[i,j]=lat[ind]
            lonb[i,j]=lon[ind]
    ga=calc_jacob(lonb,latb,P,n)
    u1,v1=ca2cs_sp(ub,vb,latb,lonb)
 
    Dx=ga*u1*Dx/a**2
    Dy=ga*v1*Dy/a**2
    return Dx,Dy

def calc_jacob(lon,lat,pn,n):
    a=1/np.sqrt(3)
    if pn==1:
        ga=cos(lon)**3*cos(lat)**3
    elif pn==2:
        clon=lon-pi/2
        ga=cos(clon)**3*cos(lat)**3
    elif pn==3:
        clon=lon.copy()
        for i in range(lon.shape[0]):
            for j in range(n):
                if clon[i,j]<0:
                    clon[i,j]=lon[i,j]+pi*2
        clon-=pi
        ga=cos(clon)**3*cos(lat)**3
    elif pn==4:
        clon=lon+pi/2
        ga=cos(clon)**3*cos(lat)**3
    elif pn==5:
        clat=lat.copy()
        ga=sin(lat)**3
    elif pn==6:
        clat=lat.copy()
        ga=sin(lat)**3
    return ga
if __name__ == '__main__':
    print('work')