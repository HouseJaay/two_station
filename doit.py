import sys
sys.path.insert(0,'/home/hao_shijie/work/two_station')
import cal_station as cs
from core import Disp,Basket
import haotool as ht
import cut
import two_station as ts

dep_max = 30
dist_max = 9000
dist_min = 2000
mag_min = 5.8

evt_full = ht.read_event('evt_2004_2016')
sta = ht.read_station('station_us')
pair = ht.mk_sta_pairs(sta)
pair.loc[:,'event'] = pair.apply(cs.do_check,axis='columns',
args=(evt_full,dep_max,dist_min,dist_max,mag_min))
pair_temp = pair[pair['event']>5]

basket = Basket(pair_temp,evt_full)

for i in basket.data.index:
    basket.data[i].prepfile('/home/hao_shijie/data/thesis/')
    cut.do_cut(basket.data[i])
    basket.data[i].filt_snr(5)
