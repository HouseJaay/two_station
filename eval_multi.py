#!/usr/bin/eny python
# evaluate events that may be interferenced by another event
# only calculate first station of station pairs
#TODO two station 
# 3.20 add pick by epicenter and depth

import pandas
from distaz import distaz
from obspy.core import UTCDateTime

VMIN = 2.5
VMAX = 7
file_pair = 'pair_temp'
file_event = 'event'
in_prefix = 'temp/'

def cal_window(stla,stlo,evla,evlo,evtime,vmin,vmax):
    dist = distaz(stla,stlo,evla,evlo).degreesToKilometers()
    start = evtime + dist/vmax
    end = evtime + dist/vmin
    return (start,end)

def read_time(evt,n):
    (year,month,day,hour,minute,sec) = (evt['year'][n],evt['month'][n],evt['day'][n],evt['hour'][n],evt['min'][n],evt['sec'][n])
    return UTCDateTime(year,month,day,hour,minute,sec)

def check_window(evt,start_p,end_p,stla,stlo,vmin,vmax):
    for j in evt.index:
        evla,evlo = evt['lat'][j],evt['lon'][j]
        time = read_time(evt,j)
        (start,end) = cal_window(stla,stlo,evla,evlo,time,vmin,vmax)
        if( abs(start-start_p) < 0.1 and abs(end-end_p) < 0.1):
            continue
        elif(not (end_p < start or start_p > end)):
            return True
    return False

def check_epdist(evla,evlo,stla,stlo):
    dist = distaz(stla,stlo,evla,evlo).degreesToKilometers()
    if(dist > 30000):
        return True
    else:
        return False

def check_dep(evdp):
    if(evdp > 100):
        return True
    else:
        return False

evt = pandas.read_table(file_event,sep='\s+')
pairs = pandas.read_table(file_pair,sep='\s+')

for np in range(len(pairs)):
    stla,stlo = pairs['lat1'][np],pairs['lon1'][np]
    in_name = in_prefix + pairs['sta1'][np] + '_' + pairs['sta2'][np] + '.lst'
    evt_p = pandas.read_table(in_name,sep='\s+')
    yes_list = []
    no_list = []
    for i in range(len(evt_p)):
        evla_p,evlo_p = evt_p['lat'][i],evt_p['lon'][i]
        time_p = read_time(evt_p,i)
        (start_p,end_p) = cal_window(stla,stlo,evla_p,evlo_p,time_p,VMIN,VMAX)
        mask = ( abs(evt['jday'] - evt_p['jday'][i]) < 2 ) & ( evt['year'] == evt_p['year'][i] )
        evt_temp = evt[mask]
        b1 = check_window(evt_temp,start_p,end_p,stla,stlo,VMIN,VMAX)
        b2 = check_epdist(evla_p,evlo_p,stla,stlo)
        b3 = check_dep(evt_p['dep'][i])
        if(b1 or b2 or b3):
            print(b1,b2,b3,i)
            no_list.append(i)
        else:
            yes_list.append(i)
    print(in_name,len(no_list),len(yes_list))
    evt_p.ix[yes_list].to_csv(in_name,sep=' ',index=False)
