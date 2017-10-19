import numpy as np
import pandas as pd
from pandas import DataFrame,Series
import os
from os.path import join
import obspy
import matplotlib.pyplot as plt
import cal_station as cs

def cal_std_mean(file_list):
    temp = np.loadtxt(file_list[0])
    temp = temp.transpose()
    t = temp[0,temp[0].argsort()]
    v = temp[1,temp[0].argsort()]
    for i in range(1,len(file_list)):
        temp = np.loadtxt(file_list[i])
        temp = temp.transpose()
        tempv = temp[1,temp[0].argsort()]
        v = np.vstack((v,tempv))
    mean_value = np.mean(v,axis=0)
    std_value = np.std(v,axis=0)
    return t,std_value,mean_value

def plot_disp_ll(file_list2,forward = '',pmin=20,pmax=60):
    fig,ax = plt.subplots()
    ax.set_xlabel('period(s)')
    ax.set_ylabel('phase velocity')
    ax.set_ylim(3,4.5)
    ax.set_xlim(pmin,pmax)
    handle = []
    for file_list in file_list2:
        t,std_v,mean_v = cal_std_mean(file_list)
        filename = os.path.basename(file_list[0])
        label = filename.split('.')[0]
        temp, = ax.plot(t,mean_v,'-',label=label)
        handle.append(temp)
    if forward:
        disp2 = np.loadtxt(forward)
        disp2 = disp2.transpose()
        temp, = ax.plot(disp2[0],disp2[1],'-',label='forward',color='black')
        handle.append(temp)
    ax.legend(handles=handle)
    plt.show()


class Disp(object):
    """
    >>> Disp(pair_row,evt)
    """
    def __init__(self, pair_row, evt):
        self.pair = pair_row.copy()
        self.evt = evt.copy()
    
    def prepfile(self,rtdir):
        CH = 'BHZ'
        def path_to_data(evt_row,rtdir,pair_row,code):
            """
            code: integer 1 or 2
            corresponding to station1 and station2 in pair
            """
            dir1 = pair_row['station1']+'_'+pair_row['station2']
            dir2 = evt_row['time'].strftime("%Y_%m_%d_%H_%M")
            temp = 'station' + str(code)
            filename = pair_row[temp] + '.' + CH
            return join(rtdir,'out',dir1,dir2,filename)

        def path_to_disp(evt_row,rtdir,pair_row):
            part1 = pair_row['station1']+'_'+pair_row['station2']
            part2 = evt_row['time'].strftime("%Y_%m_%d_%H_%M")
            dispfile = part1 + '.' + part2 + '.disp'
            return join(rtdir,'out',dispfile)
        
        self.evt.loc[:,'data1'] = self.evt.apply(path_to_data,axis='columns',args=
        (rtdir,self.pair,1))
        self.evt.loc[:,'data2'] = self.evt.apply(path_to_data,axis='columns',args=
        (rtdir,self.pair,2))

        self.evt.loc[:,'disp'] = self.evt.apply(path_to_disp,axis='columns',args=
        (rtdir,self.pair))

        path_to_sac = join(rtdir,self.pair['station1']+'_'+self.pair['station2'])
        respname1 = 'RESP.' + self.pair['net1'] + '.' + self.pair['station1']+'..'+CH
        respname2 = 'RESP.' + self.pair['net2'] + '.' + self.pair['station2']+'..'+CH
        self.evt.loc[:,'resp1'] = join(path_to_sac,respname1)
        self.evt.loc[:,'resp2'] = join(path_to_sac,respname2)
        
        self.pair.loc['raw_data'] = path_to_sac

    def check_file(self):
        """
        check if dispfiles exist
        if not, delete this event
        """
        no_file = []
        for index,e in self.evt.iterrows():
            if os.path.exists(e['disp']):
                pass
            else:
                no_file.append(index)
        self.evt = self.evt.drop(no_file)
    
    def filt_snr(self,threshold):
        def cal_snr(filename):
            st = obspy.read(filename)
            data = abs(st[0].data)
            data_max = data.max()
            data_mean = data.mean()
            return data_max/data_mean
        waste = []
        for index,e in self.evt.iterrows():
            snr1 = cal_snr(e['data1'])
            snr2 = cal_snr(e['data2'])
            if snr1>threshold and snr2>threshold:
                pass
            else:
                waste.append(index)
        print("drop %d events" % len(waste))
        self.evt = self.evt.drop(waste)

    def plot_disp(self,forward = '',pmin=20,pmax=60):
        fig,ax = plt.subplots()
        ax.set_xlabel('period(s)')
        ax.set_ylabel('phase velocity')
        ax.set_ylim(3,4.5)
        ax.set_xlim(pmin,pmax)
        handle = []
        for index,e in self.evt.iterrows():
            label = index
            disp = np.loadtxt(e['disp'])
            disp = disp.transpose()
            indexer = disp[0].argsort()
            temp, = ax.plot(disp[0][indexer],disp[1][indexer],'-',label=label)
            handle.append(temp)
        if forward:
            disp2 = np.loadtxt(forward)
            disp2 = disp2.transpose()
            temp, = ax.plot(disp2[0],disp2[1],'-',label='forward',color='black')
            handle.append(temp)

        ax.legend(handles=handle)
        plt.show()
    
    def plot_waveform(self,evt_id):
        """
        plot waveform
        evt_id is index of event in self.evt dataframe
        """
        st = obspy.read(self.evt.loc[evt_id,'data1'])
        st += obspy.read(self.evt.loc[evt_id,'data2'])
        st.plot()

    def plot_evst(self):
        pass


class Basket(object):
    def __init__(self,pairs,evt_full,evt_rule=(30,2000,9000,5.8)):
        self.pairs = pairs.copy()
        self.data = pd.Series(0,index=self.pairs.index)
        for i in self.data.index:
            self.data[i] = Disp(self.pairs.loc[i],
            cs.get_event(self.pairs.loc[i],evt_full,evt_rule[0],
            evt_rule[1],evt_rule[2],evt_rule[3]))

    def plot_mean(self,index,forward='',pmin=20,pmax=60):
        file_list2 = []
        for i in index:
            file_list2.append( list(self.data[i].evt['disp']) )

        plot_disp_ll(file_list2,forward=forward,pmin=pmin,pmax=pmax)
        

        


if __name__=='__main__':
    import doctest
    #doctest.testmod()
