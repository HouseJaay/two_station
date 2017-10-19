#!/usr/bin/env python
# TODO pick two files
from obspy.core import UTCDateTime
import obspy
import pandas
from distaz import distaz
import os
from os.path import join
from subprocess import Popen,PIPE
import sys
from glob import glob

VMAX = 7
VMIN = 2.5
freq_lim = (0.005, 0.006, 1, 1.2)
decim = (2, 5)

def pick(start,end,file_list,file):
    for filename in file_list:
        st = obspy.read(filename,headonly=True)
        b,e = st[0].stats.starttime,st[0].stats.endtime
        if( start > b and end < e ):
            file.append(filename)
            return True
    return False

def check(stats,file1,file2):
    st = obspy.read(file1,headonly=True)
    st += obspy.read(file2,headonly=True)
    if(st[0].stats.delta == st[1].stats.delta):
        stats.append(st[0].stats.delta)
        return True
    else:
        return False

def trans(file,file_resp,f,d):
    p = Popen(['sac'], stdin=PIPE, stdout=PIPE)
    s = ""
    s += "r %s\n" % file
    s += "decimate %d;decimate %d\n" % (d[0],d[1])
    s += "rmean;rtrend\n"
    s += "transfer from evalresp fname %s to none freq %f %f %f %f\n" % (file_resp,f[0],f[1],f[2],f[3])
    s += "w append .temp\n"
    s += "q\n"
    p.communicate(s.encode())

def cut(filename,start,end,dist,outfile):
    st = obspy.read(filename)
    delta = st[0].stats.delta
    width = end - start
    npts = int(width/delta)
    n = int((start - st[0].stats.starttime)/delta)
    st[0].data = st[0].data[n:n+npts]
    st[0].stats.starttime = start
    st[0].stats.sac.unused23 = int(dist)
    st[0].write(outfile,format='SAC')
    

def do_cut(Disp):
    CH = 'BHZ'
    raw_path = Disp.pair['raw_data']
    file_list1 = glob(join(raw_path,"*"+Disp.pair['station1']+'*'+CH+'*'+'SAC'))
    file_list2 = glob(join(raw_path,"*"+Disp.pair['station2']+'*'+CH+'*'+'SAC'))
    no_data = []
    for index,e in Disp.evt.iterrows():
        start = e['time'] + min(e['dist'])/VMAX
        end = e['time'] + max(e['dist'])/VMIN
        file1 = []
        file2 = []
        if(not (pick(start,end,file_list1,file1) and pick(start,end,file_list2,file2)) ):
            print("pick raw data failed %s"%e['data1'],file=sys.stderr)
            no_data.append(index)
            continue
        outdir = os.path.dirname(e['data1'])
        os.system("mkdir -p "+outdir)
        trans(file1[0],e['resp1'],freq_lim,decim)
        trans(file2[0],e['resp2'],freq_lim,decim)
        cut(file1[0] + '.temp',start,end,e['dist'][0],e['data1'])
        cut(file2[0] + '.temp',start,end,e['dist'][1],e['data2'])
    Disp.evt = Disp.evt.drop(no_data)
