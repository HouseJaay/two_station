#!/usr/bin/eny python

import pandas
from distaz import distaz
from obspy.core import UTCDateTime

VMIN = 2.5
VMAX = 7
file_pair = 'sta_pairs'
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
        flag = 0
        for j in evt_temp.index:
            evla,evlo = evt_temp['lat'][j],evt_temp['lon'][j]
            time = read_time(evt_temp,j)
            (start,end) = cal_window(stla,stlo,evla,evlo,time,VMIN,VMAX)
            if( abs(start-start_p) < 0.1 and abs(end-end_p) < 0.1):
                continue
            if( not (end_p < start or start_p > end)):
                flag = 1
                break
        if( flag == 1):
            no_list.append(i)
        else:
            yes_list.append(i)
    print(in_name,len(no_list),len(yes_list))
    evt_p.ix[yes_list].to_csv(in_name,sep=' ',index=False)
