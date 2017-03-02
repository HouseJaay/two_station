#! /usr/bin/env python
# server version 
# dev1 is cluster version

import obspy
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
import distaz

delta=0.025
npts=213163
dist=7.78
umin=2
umax=8
pmin=2
pmax=80
data1='CUT_S05C.BHZ'
data2='CUT_S06C.BHZ'

len_cor=2*npts-1
t=np.arange(1,(len_cor+1)/2)*delta
u=dist/t
mask=(u>umin)*(u<umax)
u=u[mask]
p=np.arange(pmin,pmax)
COR=np.zeros((len(p),len(u)))
U,P=np.meshgrid(u,p)

def pick(cor,uini):
    for i in range(len(u)):
        if(u[i]<uini):
            j=i
            break
    if(cor[j+1]>cor[j]):
        while(cor[j+1]>cor[j]):
            j+=1
        i=j
    elif(cor[j-1]>cor[j]):
        while(cor[j-1]>cor[j]):
            j-=1
        i=j
    return u[i]    

def onclick(event):
    global click_x,click_y
    click_x,click_y=event.xdata,event.ydata


st=obspy.read(data1)  # attention: sequence
st+=obspy.read(data2) 
row=0
for period in range(pmin,pmax):
    b=signal.firwin(501,[1.0/(period+0.2),1.0/(period-0.2)],window=('kaiser',9),nyq=1/delta/2,pass_zero=False)
# cut window
   #st[0].data=st[0].data*w1
    #st[1].data=st[1].data*w2
# filter
    array1=signal.lfilter(b,1,st[0].data)
    array2=signal.lfilter(b,1,st[1].data)
# normalize
    array2=array2*array1.max()/array2.max()
    #st[0].write('g11_bp20uw',format='SAC')
    #st[1].write('g10_bp20uw',format='SAC')
# correlate , first input signal has larger epicenter distance
    corr=signal.correlate(array1,array2,mode='full')
# data prepare
    cor=corr[(len_cor+1)/2:len_cor]
    cor=cor[mask]
    cor=cor/cor.max() #normalize
    print(cor.max())
    COR[row]=cor
    row+=1
# plot correlation function   
    #plt.figure()
    #plt.plot(u[20:100],cor[20:100],'b')
    #plt.show()
# pick
    #i=pick(u,cor,Uini)
    #Uini=u[i]
    #print(str(period)+' '+str(u[i]))

fig,ax=plt.subplots()
cf=ax.contourf(P,U,COR)
fig.colorbar(cf)
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

pini=int(click_x)
uini=pick(COR[pini-pmin],click_y)
print(pini,uini)
utemp=uini
for period in range(pini+1,pmax,1):
    utemp=pick(COR[period-pmin],utemp)
    print(period,utemp)
utemp=uini
for period in range(pini-1,pmin-1,-1):
    utemp=pick(COR[period-pmin],utemp)
    print(period,utemp)
