import xarray as xr
import numpy as np

# Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):
def calc_z(typ,ds):
    # Function to calculate depth of levels from CROCO output. This script has been adapted
    # from the croco_tools/Preprocessing_tools/zlevs.m script.
    # 
    # Inputs:
    # 
    # ds = CROCO history or average file xarray dataset
    # typ = 'r'; rho-points, 'w'; w-points
    #
    # Outputs:
    # z = depth of rho or w points
    
    if(ds.Vtransform != 2):
        print('ERROR - wrong Vtransform')
        return
    
    if (typ=='r'):
        Cs = ds.Cs_r
        sc = ds.sc_r
        z = xr.zeros_like(ds.temp).rename('z_rho')
        N = len(ds.s_rho)
    elif (typ=='w'):
        Cs = ds.Cs_w
        sc = ds.sc_w
        z = xr.zeros_like(ds.w).rename('z_w')
        N = len(ds.s_w)
    
    h = ds.h#.where(ds.h==0.,1.e-2,ds.h)
    #Dcrit = 0.01
    zeta = ds.zeta#.where(ds.zeta<(Dcrit-h),Dcrit-h,ds.zeta)
    
    hinv=1/h;
    h2=(h+ds.hc)
    cff = ds.hc*sc
    h2inv = 1/h2
    
    z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
    return(z)

# Define a function to calculate depths (python script adapted from croco_tools/Preprocessing_tools/zlevs.m):
def calc_z1D(typ,ds,h,zeta):
    # Function to calculate 1D depth of levels from CROCO output given an input H and zeta
    # 
    # Inputs:
    # 
    # ds = CROCO history or average file xarray dataset
    # typ = 'r'; rho-points, 'w'; w-points
    #
    # Outputs:
    # z = depth of rho or w points
    
    if(ds.Vtransform != 2):
        print('ERROR - wrong Vtransform')
        return
    
    if (typ=='r'):
        Cs = ds.Cs_r
        sc = ds.sc_r
        z = xr.zeros_like(Cs).rename('z_rho')
        N = len(ds.s_rho)
    elif (typ=='w'):
        Cs = ds.Cs_w
        sc = ds.sc_w
        z = xr.zeros_like(Cs).rename('z_w')
        N = len(ds.s_w)
    
    hinv=1/h;
    h2=(h+ds.hc)
    cff = ds.hc*sc
    h2inv = 1/h2
    
    z = (cff+Cs*h)*h/h2 + zeta*(1+(cff+Cs*h)*h2inv)
    
    return(z)


