#!/usr/bin/env python

import obspy
from scipy import signal
import numpy as np

file = ('S05C.BHZ','S06C.BHZ')
flag = True
prange = (20,150)

def cal_snr(file,w_sac,pr):
    st = obspy.read(file[0])
    st += obspy.read(file[1])
    delta = st[0].stats.delta
    f1 = open(file[0]+'.snr','w')
    f2 = open(file[1]+'.snr','w')
    for period in range(pr[0],pr[1],10):
        b = signal.firwin(1001,[1.0/(period+0.2),1.0/(period-0.2)],window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
        cp = st.copy()
        cp[0].data = signal.lfilter(b,1,st[0].data)
        cp[1].data = signal.lfilter(b,1,st[1].data)
        if( w_sac ):
            cp[0].write(file[0]+'_'+str(period),format='SAC')
            cp[1].write(file[1]+'_'+str(period),format='SAC')

        snr1 = cp[0].data.max()/np.absolute(cp[0].data).mean()
        snr2 = cp[1].data.max()/np.absolute(cp[1].data).mean()
        f1.write("%d %f\n" %(period,snr1))
        f2.write("%d %f\n" %(period,snr2))
    f1.close()
    f2.close()

cal_snr(file,flag,prange)
