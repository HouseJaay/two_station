#!/usr/bin/env python
# input file: event list , station pairs list
# output file: selected event lists of station pairs

import pandas
import distaz
filename_sta='pair_temp'
filename_evt='event'
out_prefix="temp/"

sta = pandas.read_table(filename_sta,sep='\s+')
evt = pandas.read_table(filename_evt,sep='\s+')
out_header=['year','month','day','jday','hour','min','sec','lat','lon','dep','mw']

for npair in range(len(sta)):
    row_list=[]
    stla1,stlo1 = sta['lat1'][npair],sta['lon1'][npair]
    stla2,stlo2 = sta['lat2'][npair],sta['lon2'][npair]
    sta_az = distaz.distaz(stla2,stlo2,stla1,stlo1).getAz()
    name_pairs = out_prefix+sta['sta1'][npair]+'_'+sta['sta2'][npair]+'.lst'
    for nevt in range(len(evt)):
        evla = evt['lat'][nevt]
        evlo = evt['lon'][nevt]
        evt_az = distaz.distaz(evla,evlo,stla1,stlo1).getAz()
        if(abs(sta_az-evt_az)<2.0 or abs(abs(sta_az-evt_az)-180)<2.0):
            row_list.append(nevt)
    evt.ix[row_list].to_csv(name_pairs,sep=' ',header=out_header,index=False)
