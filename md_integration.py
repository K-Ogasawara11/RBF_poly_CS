import numpy as np
from scipy.spatial import KDTree
from numpy import pi,sin,cos,tan
from  params import grav,rea

def calc_dh_adv(SDx,SDy,u,v,p_vec2all,h):
    pdh=(SDx@h)*u
    pdh+=(SDy@h)*v
    dh=p_vec2all@pdh
    return dh


def calc_tendencyv2(Dx,Dy,hiv,cui,cvi,ui,vi,hi,rgi,col,dt):
    #calc_vorticity
    vor=((Dx@cvi)-(Dy@cui))/rgi
    
    #calc static energy
    E=(cui*ui+cvi*vi)/2+grav*hi
   

    dh=-(Dx@(rgi*ui*hi)+Dy@(rgi*vi*hi))/rgi+hiv@hi   
    du=-Dx@E+rgi*vi*(col+vor)+hiv@cui
    dv=-Dy@E-rgi*ui*(col+vor)+hiv@cvi
    return du*dt,dv*dt,dh*dt

def calc_tendency_surf(Dx,Dy,hiv,cui,cvi,ui,vi,u,v,hi,hs,rg,rgi,col,dt):
    #calc geopotential
    gp=hi+hs

    #calc_vorticity
    vor=((Dx@cvi)-(Dy@cui))/rgi
    
    #calc static energy
    E=(cui*ui+cvi*vi)/2+grav*gp
   

    dh=-(Dx@(rgi*ui*hi)+Dy@(rgi*vi*hi))/rgi+hiv@hi
   
    du=-Dx@E+rgi*vi*(col+vor)+hiv@cui
    dv=-Dy@E-rgi*ui*(col+vor)+hiv@cvi
    return du*dt,dv*dt,dh*dt


def calc_tendency_adv(Dx,Dy,hiv,u,v,h,dt):
    dh=-u*(Dx@h)-v*(Dy@h)+hiv@h
    return dh*dt


if __name__ == '__main__':
    print('work')