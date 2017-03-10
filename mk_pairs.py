#!/usr/bin/env python
# read list of names and output pairs of names

filename='sta'
outname='sta_pairs'

with open(filename,'r') as f:
    mylist = f.read().splitlines()

outlist=[]
for i in range(len(mylist)-1):
    for j in range(i+1,len(mylist)):
        outlist.append(mylist[i]+' '+mylist[j]+'\n')

with open(outname,'w') as f:
    f.write("sta1 lat1 lon1 sta2 lat2 lon2\n"+"".join(outlist))
