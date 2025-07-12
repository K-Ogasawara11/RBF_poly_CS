import numpy as np
from numba import jit
from scipy.spatial import KDTree
from numpy import pi,sin,cos,tan
import scipy
from params import rea

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
    #r=np.zeros((sx.size,n))
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


def calc_Dmat_poly_eang(RA,R,r,dx,dy,dz,stx,sty,bdx,bdy,sx,n,P):
    a=1/np.sqrt(3)
    Dx=np.zeros((sx.size,n))
    Dy=np.zeros((sx.size,n))
    alpha=np.arctan(stx/a)
    beta=np.arctan(sty/a)
    xx,xy,xz,yx,yy,yz=calc_coefficient_eang(a,stx,sty,r,alpha,beta,P)
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

def calc_coefficient_eang(a,x,y,r,alpha,beta,pn):
    r3=r**3;r2=r**2
    x2=x**2;y2=y**2
    if pn==1:
        xx=-a**2*x/(r3*(cos(alpha)**2));xy=a*(r2*cos(alpha)**2-x2)/(r3*(cos(alpha)**2));xz=-y*x*a/(r3*(cos(alpha)**2))
        yx=-a**2*y/(r3*cos(beta)**2);yy=-y*x*a/(r3*cos(beta)**2);yz=a*(r2*cos(beta)**2-y2)/(r3*cos(beta)**2)
      
    elif pn==2:
        xx=a*(-r2*cos(alpha)**2+x2)/(r3*(cos(alpha)**2));xy=-a**2*x/(r3*(cos(alpha)**2));xz=-y*x*a/(r3*(cos(alpha)**2))
        yx=x*y*a/(r3*cos(beta)**2);yy=-a**2*y/(r3*cos(beta)**2);yz=a*(r2*cos(beta)**2-y2)/(r3*cos(beta)**2)
        
       
    elif pn==3:
        xx=a**2*x/(r3*(cos(alpha)**2));xy=a*(-r2*cos(alpha)**2+x2)/(r3*(cos(alpha)**2));xz=-a*x*y/(r3*(cos(alpha)**2))
        yx=a**2*y/(r3*cos(beta)**2);yy=(a*y*x)/(r3*cos(beta)**2);yz=a*(r2*cos(beta)**2-y2)/(r3*cos(beta)**2)
   
    elif pn==4:
        xx=a*(r2*cos(alpha)**2-x2)/(r3*(cos(alpha)**2));xy=a**2*x/(r3*(cos(alpha)**2));xz=-a*y*x/(r3*(cos(alpha)**2))
        yx=-x*y*a/(r3*cos(beta)**2);yy=a**2*y/(r3*cos(beta)**2);yz=a*(r2*cos(beta)**2-y2)/(r3*cos(beta)**2)

    elif pn==5:
        xx=x*y*a/(r3*(cos(alpha)**2));xy=a*(r2*cos(alpha)**2-x2)/(r3*(cos(alpha)**2));xz=-a**2*x/(r3*(cos(alpha)**2))
        yx=a*(-r2*cos(beta)**2+y2)/(r3*cos(beta)**2);yy=-x*y*a/(r3*cos(beta)**2);yz=-a**2*y/(r3*cos(beta)**2)
    elif pn==6:
        xx=(-x*y*a)/(r3*(cos(alpha)**2));xy=a*(r2*cos(alpha)**2-x2)/(r3*(cos(alpha)**2));xz=a**2*x/(r3*(cos(alpha)**2))
        yx=a*(r2*cos(beta)**2-y2)/(r3*cos(beta)**2);yy=-x*y*a/(r3*cos(beta)**2);yz=a**2*y/(r3*cos(beta)**2)

    return xx,xy,xz,yx,yy,yz

def calc_Dmatrix_phs(x,y,z,x1,y1,z1,idim,pn,n,deg,projection):
    index=np.zeros((idim.size,n))
    tree = KDTree(np.array([x1,y1,z1]).T)
    for i in range(idim.size):
        ind=int(idim[i])
        value,indexb=tree.query([x[ind],y[ind],z[ind]],k=n)
        index[i,:]=np.array(indexb)
    #ここに手を加える必要がある。
    #indexは微分行列とかけるベクトルに対する引数
    R,dx,dy,dz,stx,sty,stx0,sty0,r=calc_R_phs(idim,x1,y1,z1,index,n,pn)
    A,RA,dxpols,dypols=calc_ARA_phs(R,idim.size,n,stx0,sty0,deg,projection) 
    if projection=='edist':
        Dx,Dy=calc_Dmat_poly(RA,R,r,dx,dy,dz,stx,sty,dxpols,dypols,idim,n,pn)
    elif projection=='eang':
        Dx,Dy=calc_Dmat_poly_eang(RA,R,r,dx,dy,dz,stx,sty,dxpols,dypols,idim,n,pn)
    

    return Dx,Dy,r,index,R
    

def convert_sparse_c2inner(index,Dx,Dy,idim,N,n,vsize,c2innerL):
    row,col,vals=make_row_col_c2input(index,Dx,n,idim.size,c2innerL)
    SDx=scipy.sparse.csr_matrix((vals,(row,col)),shape = (vsize,vsize))
    row,col,vals=make_row_col_c2input(index,Dy,n,idim.size,c2innerL)
    SDy=scipy.sparse.csr_matrix((vals,(row,col)),shape = (vsize,vsize))
    return SDx,SDy


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

@jit(nopython=True)
def make_row_col_c2input(index,vals,n,N,c2innerL):
    val2=vals.reshape(N*n)
    col=index.reshape(N*n)
    count=0
    row=np.zeros(N*n)
    for i in range(N):
        for j in range(n):
            row[count]=c2innerL[i]
            count+=1
    return row,col,val2


def laplasian_fd2d(k,ep,R,N,n):
    import math
    import numpy as np
    import scipy.special as sp
    laplas2=np.zeros((N,n))
    kk=math.factorial(k)#kの階乗
    ep2=ep**2
    P=sp.genlaguerre(k,0)
    LA=np.zeros((n+1,n+1))
    BB=np.zeros(n+1)
    LA[n,:]=1
    LA[:,n]=1
    LA[n,n]=0
    for i in range(N):
        ra=R[i,:,:]
        rb=R[i,0,:]
        p=(-4)**k*kk*P(ep2*rb)
        laplas=ep2**k*p
        BB[0:n]=(laplas*np.exp(-ep2*rb)).T
        LA[0:n,0:n]=np.exp(-ep2*ra)
        lap=np.linalg.solve(LA,BB).T#;print(lap.shape)
        laplas2[i,:]=lap[0:n]
    return laplas2


def calc_d2d2(N,n,ep,index,R):
    import numpy as np
    
    ep4=ep**4
    ep2=ep**2
    A=np.exp(-ep**2*R)
    RA=np.zeros(((N,n+1,n+1)))
    RA[:,0:n,0:n]=A
    xd=np.zeros((N,n));yd=np.zeros((N,n));zd=np.zeros((N,n))
    xtx=np.zeros(((N,n,n)))
    Dx=np.zeros((N,n+1));DxT=np.zeros(n+1);Dxb=np.zeros(n+1)
    #Dxb=np.zeros((n+1,n+1));Dyb=np.zeros((n+1,n+1));Dzb=np.zeros((n+1,n+1))
    for i in range(N):
        RA[i,n,0:n]=1;RA[i,0:n,n]=1
        RA[i,n,n]=0
        Dxb[0:n] = (4*ep4*(R[i,0,:])-4*ep2)*RA[i,0,0:n] #;print(Dxb.shape)
        Dx[i,:]=np.linalg.solve(RA[i,:,:],Dxb)
    return Dx[:,0:n]/rea**2




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


def spherical_laplas(r,ep,idim,n):
    ep2=ep**2
    ep4=ep2**2
    lap=np.zeros((idim.size,n))
    RA=np.zeros((n+1,n+1))
    RA[n,:]=1;RA[:,n]=1;RA[n,n]=0
    B=np.zeros(n+1)
    for i in range(idim.size):
        rb=r[i,0,:]
        RA[0:n,0:n]=np.exp(-ep2*r[i,:,:])
        phi=np.exp(-ep2*rb)
        phir=-2*ep2*phi
        phir2=(4*ep4*rb-2*ep2)*phi
        lapb1=0.25*((4-rb)*phir2+(4-3*rb)*phir)
        B[0:n]=lapb1
        lapb=np.linalg.solve(RA,B)
        lap[i,:]=lapb[0:n]
    return lap



if __name__ == '__main__':
    print('work')