#!/usr/bin/env python

import haotool as ht
from distaz import distaz
import numpy as np

# mean distance of station pairs 
def cal_dist(row_evt,row_pair):
    evla,evlo = row_evt['lat'],row_evt['lon']
    stla1,stlo1 = row_pair['lat1'],row_pair['lon1']
    stla2,stlo2 = row_pair['lat2'],row_pair['lon2']
    dist1 = distaz(stla1,stlo1,evla,evlo).degreesToKilometers()
    dist2 = distaz(stla2,stlo2,evla,evlo).degreesToKilometers()
    return dist1,dist2

def check_az(row_evt,row_pair,sta_az):
    evla,evlo = row_evt['lat'],row_evt['lon']
    stla1,stlo1 = row_pair['lat1'],row_pair['lon1']
    evt_az = distaz(evla,evlo,stla1,stlo1).getAz()
    if(abs(sta_az-evt_az)<2.0 or abs(abs(sta_az-evt_az)-180)<2.0):
        return True
    else:
        return False

def cal_window(row_evt,row_pair,vmin,vmax):
    dist1,dist2 = cal_dist(row_evt,row_pair)
    evtime = row_evt['time']
    start1 = evtime + dist1/vmax
    end1 = evtime + dist1/vmin
    start2 = evtime + dist2/vmax
    end2 = evtime + dist2/vmin
    return start1,end1,start2,end2

def check_multi(row_evt,row_pair,evt_full):
    vmin,vmax = 2.5,7
    t1 = cal_window(row_evt,row_pair,vmin,vmax)
    mask = (abs(evt_full['time'] - row_evt['time']) < 86400)
    evt_full = evt_full[mask]
    t2 = evt_full.apply(cal_window,axis='columns',args=(row_pair,vmin,vmax))
    def check_window(start1,end1,start2,end2):
        if(abs(start1-start2) < 0.1 and abs(end1-end2) < 0.1):
            return True
        elif(not (start1 > end2 or start2 > end1)):
            return False
        else:
            return True
    for i in t2.index:
        b1 = check_window(t1[0],t1[1],t2[i][0],t2[i][1])
        b2 = check_window(t1[2],t1[3],t2[i][2],t2[i][3])
        if(not (b1 and b2)):
            return False
    return True


def do_check(row_pair,evt,dep_max,dist_min,dist_max,mag_min):

    start,end = row_pair['start'],row_pair['end']
    mask0 = (evt['time'] > start) & (evt['time'] < end - 86400)
    evt1 = evt[mask0].copy()

    evt2 = evt1[evt1['dep'] <= dep_max].copy()
    evt2 = evt2[evt2['mw'] > mag_min]
    
    stla1,stlo1 = row_pair['lat1'],row_pair['lon1']
    stla2,stlo2 = row_pair['lat2'],row_pair['lon2']
    sta_az = distaz(stla2,stlo2,stla1,stlo1).getAz()
    mask1 = evt2.apply(check_az,axis='columns',args=(row_pair,sta_az))
    evt2 = evt2[mask1]

    if(len(evt2) == 0):
        return 0
    evt2.loc[:,'dist'] = evt2.apply(cal_dist,axis='columns',args=(row_pair,))
    mean = lambda x: (x[0]+x[1])/2.0
    evt2 = evt2[ (evt2['dist'].map(max) < dist_max) & (evt2['dist'].map(min) > dist_min)]

    mask2 = evt2.apply(check_multi,axis='columns',args=(row_pair,evt1))
    evt2 = evt2[mask2]
    
    return len(evt2)

def get_event(row_pair,evt,dep_max,dist_min,dist_max,mag_min):

    start,end = row_pair['start'],row_pair['end']
    mask0 = (evt['time'] > start) & (evt['time'] < end - 86400)
    evt1 = evt[mask0].copy()

    evt2 = evt1[evt1['dep'] <= dep_max].copy()
    evt2 = evt2[evt2['mw'] > mag_min]
    stla1,stlo1 = row_pair['lat1'],row_pair['lon1']
    stla2,stlo2 = row_pair['lat2'],row_pair['lon2']
    sta_az = distaz(stla2,stlo2,stla1,stlo1).getAz()
    mask1 = evt2.apply(check_az,axis='columns',args=(row_pair,sta_az))
    evt2 = evt2[mask1]

    if(len(evt2) == 0):
        return 0
    evt2.loc[:,'dist'] = evt2.apply(cal_dist,axis='columns',args=(row_pair,))
    mean = lambda x: (x[0]+x[1])/2.0
    evt2 = evt2[ (evt2['dist'].map(max) < dist_max) & (evt2['dist'].map(min) > dist_min)]

    mask2 = evt2.apply(check_multi,axis='columns',args=(row_pair,evt1))
    evt2 = evt2[mask2]
    
    return evt2

