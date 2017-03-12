#! /usr/bin/env python

import obspy
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import distaz
import pandas
from obspy.core import UTCDateTime

sta_pairs = 'pair_temp'
VRANGE = (2,8)
PRANGE = (20,100)
CH = 'BHZ'

def read_time(evt,n):
    (year,month,day,hour,minute,sec) = (evt['year'][n],evt['month'][n],evt['day'][n],evt['hour'][n],evt['min'][n],evt['sec'][n])
    return UTCDateTime(year,month,day,hour,minute,sec)

def read_pairs(pair,n,ch):
    stla1,stlo1 = pair['lat1'][n],pair['lon1'][n]
    stla2,stlo2 = pair['lat2'][n],pair['lon2'][n]
    dist = distaz.distaz(stla1,stlo1,stla2,stlo2).degreesToKilometers()
    sac_path = 'out/' + pair['sta1'][n] + '_' + pair['sta2'][n] + '/'
    evt_file = 'temp/' + pair['sta1'][n] + '_' + pair['sta2'][n] + '.lst'
    evt = pandas.read_table(evt_file,sep='\s+')
    filelist = []
    fileout = []
    for i in range(len(evt)):
        evtime = read_time(evt,i)
        evt_dir = evtime.strftime("%Y_%m_%d_%H") + '/'
        file1 = sac_path + evt_dir + pair['sta1'][n] + '.' + ch
        file2 = sac_path + evt_dir + pair['sta2'][n] + '.' + ch
        filelist.append((file1,file2))
        fileout.append('out/' + pair['sta1'][n] + '_' + pair['sta2'][n] + '.' + evtime.strftime("%Y_%m_%d_%H") + '.disp')
    return dist,filelist,fileout

def pick(cor,uini,u):
    for i in range(len(u)):
        if(u[i]<=uini):
            j=i
            break
    if(j==0 or j==(len(u)-1)):
        return -1
    if(cor[j+1]>cor[j]):
        while(j<(len(u)-1) and cor[j+1]>cor[j]):
            j+=1
        i=j
    elif(cor[j-1]>cor[j]):
        while(j>0 and cor[j-1]>cor[j]):
            j-=1
        i=j
    return u[i]    

def onclick(event):
    global click_x,click_y
    click_x,click_y=event.xdata,event.ydata


#file2 = 'out/S05C_S06C/2006_02_04_14/S05C.BHZ'
#file1 = 'out/S05C_S06C/2006_02_04_14/S06C.BHZ'
#file_out = 'out/S05C_S06C.2006_02_04_14.disp'
def two_station(file1,file2,dist,vrange,prange,file_out):
    st = obspy.read(file1)  # attention: sequence large distance first
    st += obspy.read(file2) 
    delta = st[0].stats.delta
    npts = st[0].stats.npts

    len_cor = 2*npts-1
    t = np.arange(1,int((len_cor+1)/2))*delta
    v = dist/t
    mask = (v>vrange[0]) * (v<vrange[1])
    v = v[mask]
    p = np.arange(prange[0], prange[1])
    COR = np.zeros((len(p),len(v)))
    V,P = np.meshgrid(v,p)

    row=0
    for period in range(prange[0],prange[1]):
        b = signal.firwin(1001,[1.0/(period+0.2),1.0/(period-0.2)],window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
# filter
# TODO sequence
        if(st[0].stats.sac.unused23 > st[1].stats.sac.unused23):
            array1 = signal.lfilter(b,1,st[0].data)
            array2 = signal.lfilter(b,1,st[1].data)
        else:
            array2 = signal.lfilter(b,1,st[0].data)
            array1 = signal.lfilter(b,1,st[1].data)
        """
        test=obspy.read(file1)
        test+=obspy.read(file2)
        test[0].data=signal.lfilter(b,1,test[0].data)
        test[1].data=signal.lfilter(b,1,test[1].data)
        test[0].write(file1+'_'+str(period),format='SAC')
        test[1].write(file2+'_'+str(period),format='SAC')
        """
# normalize
        array1 = array1/array1.max()
        array2 = array2/array2.max()
# correlate , first input signal has larger epicenter distance
        corr = signal.correlate(array1,array2,mode='full')
# data prepare
        cor = corr[int((len_cor+1)/2):len_cor]
        cor = cor[mask]
        cor = cor/cor.max() #normalize
        COR[row] = cor
        row+=1
# pick
    fig,ax = plt.subplots()
    cf = ax.contourf(P,V,COR)
    fig.colorbar(cf)
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    f = open(file_out,'w')
    pini = int(click_x)
    uini = pick(COR[pini-prange[0]], click_y, v)
    f.write("%d %f\n" % (pini,uini))
    utemp = uini
    for period in range(pini+1, prange[1] ,1):
        utemp = pick(COR[period-prange[0]], utemp, v)
        if(utemp > 0):
            f.write("%d %f\n" % (period,utemp))
        else:
            break
    utemp=uini
    for period in range(pini-1,prange[0]-1,-1):
        utemp = pick(COR[period-prange[0]], utemp, v)
        if(utemp > 0):
            f.write("%d %f\n" % (period,utemp))
        else:
            break
    f.close()

#two_station(file1, file2, 73, VRANGE, PRANGE, file_out)


pair = pandas.read_table(sta_pairs,sep='\s+')
for n in range(len(pair)):
    (dist, fl, fo) = read_pairs(pair, n, CH)
    for ne in range(len(fl)):
        two_station(fl[ne][0], fl[ne][1], dist, VRANGE, PRANGE, fo[ne])
