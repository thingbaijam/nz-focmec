"""
Set of functions for retrieving focal mechanism from the models,
and also, for any validating analysis. 

This module implements:

def get_focmec(lon=None, lat=None, dep=None, Mw = None, regime= 'crust', \
               preferred_model = 'all', 
               isgetpdf=True, modelfolder = '../models/',
               config_file = 'nzfocmec_v1.ini',
               subduction_zone = None,)
               
               

"""

import json
import toml  
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


#-------------------------------------------------------------------------------
def get_focmec(lon=None, lat=None, dep=None, Mw = None, regime= 'crust', \
               preferred_model = 'all', 
               isgetpdf=True, modelfolder = '../models/',
               config_file = 'nzfocmec_v1.ini',
               subduction_zone = None,):
    # load config
    config = toml.load(modelfolder+config_file)
    config['modelfolder'] = modelfolder
    regime = regime.lower()
    # if-else is faster than select
    if regime in ['crust']:
        if (lon is None) | (lat is None):
            raise Exception("Lon and Lat required")
        # locate ntdomain
        dom = get_ntdomain(lon, lat, config)
        with open(modelfolder+config[regime]['file']) as f:
            focmod = json.load(f)
        # get sdrp
        if dom is None:
            raise Exception("Lon, Lat is outside the neotectonic domains.")
        else:
            sdrp = get_sdrp_crust(focmod[dom], Mw, prefmod = preferred_model)    
        return sdrp    
            
    elif regime=='slab':
        if dep is None:
            raise Exception("Intraslab event requires depth (km)")
            
        # identify the subduction zone 
        if subduction_zone is None:
            if (lon is None) | (lat is None):
                raise Exception("Lon and Lat required")
            
        subduction_zone, sub_mod = get_subduction(lon, lat, config)
        
        with open(modelfolder+config[regime]['file']) as f:
            focmod = json.load(f)
        
        slabmod = focmod[subduction_zone]
        edep = get_depbincenter(dep)
        sdrp = slabmod[edep]
        if preferred_model in ['mean_all']:
            temp = sdrp[1]
            sdrp = [[temp[0], temp[1], temp[2], 1.0]]
        return sdrp
    
    elif regime=='interface':
        # identify subduction zone 
        subduction_zone, sub_mod = get_subduction(lon, lat, config)      
        fm_strike = sub_mod[0]
        fm_dip = sub_mod[1]
        fm_loc = np.transpose([[lon], [lat]])
        s = round(fm_strike(fm_loc)[0])
        d = round(fm_dip(fm_loc)[0])
        # constrained
        d = 10. if d<10. else d
        d = 45. if d>45. else d
        
        if subduction_zone == 'puy':
            s = 40. if s>40. else s
            s =  0. if s<0. else s    
        elif subduction_zone == 'hik':  
            s = 260. if s>260 else s
            s = 240. if s<240 else s
        else:
            raise("Subduction should be either hik or puy")
        return [[s,d,90,1.0]]
    else:    
        raise Exception("Regime should be either 'crust','slab' or 'interface'")
        
#---------------------------------------------------------------------------
def get_depbincenter(kdep):
    fbs, min_x,max_x = 5, 20, 300
    x_bin = [d for d in range(min_x, max_x, 2*fbs)]
    for db in x_bin:
        if (kdep>=(db-fbs)) & (kdep<(db+fbs)):
            return(str(db))
    return None
#---------------------------------------------------------------------------
def get_subduction(lon, lat, config):
    point = Point(lon, lat)
    hik_npfile = config['modelfolder'] \
                     +config['interface']['hikurangi']['file'] 
    hik_mod = np.load(hik_npfile, allow_pickle=True)[()]
    hikbounds = hik_mod[2]
    if hikbounds.contains(point):
        return 'hik', hik_mod
    else:
        puy_npfile = config['modelfolder']\
                        +config['interface']['puysegur']['file'] 
        puy_mod = np.load(hik_npfile, allow_pickle=True)[()]
        puybounds = puy_mod[2]
        if puybounds.contains(point):
            return 'puy', puy_mod
        else:
            raise Exception("Lon,Lat not in subduction zone")
# for crustal events ------------------------------------------------------------------
def get_ntdomain(elon, elat, config):
    # for a given lon, lat, return the domain
    file = config['modelfolder'] + config['ntdomains']['file']
    with open(file) as f:
        ntdomains = json.load(f)
        
    point, domain = Point(elon, elat), None

    for dn in ntdomains.keys():
        lons, lats = ntdomains[dn]['lon'], ntdomains[dn]['lat']
        dpoints = [(x, y) for x, y in zip(lons, lats)]
        polygon = Polygon(dpoints)
        if polygon.contains(point):
            domain = dn
            break;
    return domain

#---------------------------------------------------------------------------
def get_sdrp_crust(fmodel, mag, prefmod = 'mean_all'):
    # Returns a list of strike,dip, rake and probability 
    # Here, prefmod can either one of following:
    # 'case1', 'case2', 'mean_case1', 'mean_case2', 'all', 'mean_all'
    
    def isequal_fm(fm1, fm2):
        #fm = sdrp 
        for x1, x2 in zip(fm1, fm2):
            if x1!=x2: return False
        return True
   
    if prefmod in ['case1','case2']:
        tag = '>15' if mag<=7.5 else '>45'
        if tag not in fmodel[prefmod].keys():
            tag = '>15'      
        #
        strike = fmodel[prefmod][tag]['strikeAn']
        dip = fmodel[prefmod][tag]['dipAn']
        rake = fmodel[prefmod][tag]['rakeAn']
        prob = fmodel[prefmod][tag]['prob']
        sdrp =  [[s,d,r,p] for s,d,r,p in zip(strike, dip, rake, prob)]
        return sdrp
    
    elif prefmod in ['mean_case1','mean_case2']:
        # flength_tag = ['>15', '>15', '>45']
        pmod = prefmod[-5:]  # i.e. "case1" if "mean_case1"
        tag = '>15' if mag<=7.5 else '>45'
        if tag not in fmodel[pmod].keys():
            tag = '>15'  
        #
        strike = fmodel[pmod][tag]['strikeAn']
        dip = fmodel[pmod][tag]['dipAn']
        rake = fmodel[pmod][tag]['rakeAn']
        prob = fmodel[pmod][tag]['prob']
        idx = prob.index(max(prob))
        sdrp = [[strike[idx], dip[idx], rake[idx], 1.0],]
        return sdrp
    
    elif prefmod in ['all','mean_all']:
        if prefmod=='all':
            sdrp1 = get_sdrp_crust(fmodel, mag, prefmod = 'case1')
            sdrp2 = get_sdrp_crust(fmodel, mag, prefmod = 'case2')
        else:
            sdrp1 = get_sdrp_crust(fmodel, mag, prefmod = 'mean_case1')
            sdrp2 = get_sdrp_crust(fmodel, mag, prefmod = 'mean_case2')
        
        if isequal_fm(sdrp1, sdrp2):
            sdrp = sdrp1
        else:
            sdrp = []
            for x in sdrp2:
                sdrp1.append(x)
            for x in sdrp1:
                sdrp.append([x[0],x[1], x[2], x[3]*0.5]) 
        
        return sdrp
    
    else:
        raise Exception("prefmod  should be either 'case1', \
                'case2','mean_case1','mean_case2','all', or 'mean_all'")
