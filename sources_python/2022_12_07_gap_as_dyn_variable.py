
#########################################################
# Differences with my calibTraj code:
# (1)
# I calculate leading speeds from vData and gapData
# here, directly column vLead taken. Since vLead not always
# local or platoon consistent, deviations
# (strong for d8_0900_0930_road2_veh145.FCdata, much higher SSE here,
# minor in most other cases
# (2)
# had imposed lower bound GAP_MIN=0.4 m; now also here
#########################################################

import numpy as np
import pandas as pd
import math
from scipy.optimize import minimize, rosen
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
#from tqdm import tqdm
#from joblib import Parallel, delayed
#import multiprocessing

GAP_MIN=0.4  # pseudo-constant; value as in my c++ program

class IDM:
    def __init__(self,v0,T,s0,a,b):
        self.v0=v0
        self.T=T
        self.s0=s0
        self.a=a
        self.b=b

        self.speed_limit = 1000
        self.bmax=9
    
    #free acceleration equation
    '''
    @param v: actual speed (m/s)
    @return : free acceleration (m/s**2)
    '''

    
    def calcAccFree(self,v):

        # determine valid local v0

        v0eff = np.maximum(0.01,np.minimum(self.v0,self.speed_limit))

        accFree = self.a*(1-math.pow(v/v0eff,4)) if v<v0eff else self.a*(1-(v/v0eff))

        return  accFree

    # interaction Acceleration equation

    '''
    @param s:     actual gap [m]
    @param v:     actual speed [m/s]
    @param vl:    leading speed [m/s]
    @return:  acceleration [m/s^2]
    '''

    def calcAccInt(self,s,v,vl):

        sstar = self.s0 + np.maximum(0,v*self.T + 0.5*v*(v-vl)/np.sqrt(self.a*self.b))

        accInt = -self.a*math.pow(sstar/np.maximum(s,0.1*self.s0),2)

        #return np.maximum(-self.bmax,accInt)
        return accInt

    
    # Final longitudinal acceleration equation

    def calcAccLong(self,s,v,vl):
        accLong=np.maximum(-self.bmax,self.calcAccFree(v)+self.calcAccInt(s,v,vl))
        #if self.v0>6.8 and self.v0<7.0:
        if False:
            sstar = self.s0 + np.maximum(0,v*self.T + 0.5*v*(v-vl)/np.sqrt(self.a*self.b))
            print(f' calcAcclong: s=',s,' sstar=',sstar,' accIDM=',accLong)
            
            
        return accLong
  

    
def sim(v0,T,s0,a,b, FCdata):
    
    data=FCdata
    dt=0.2
    #print(f'simulate CF pair with v0={v0},T={T},s0={s0},a={a},b={b}')
    CF=IDM(v0,T,s0,a,b)    # IDM model as defined above
    

    # convert pandas dataframe to arrays since left-assignment faster
    
    #x1=np.empty(len(data), dtype=float)
    v1=np.empty(len(data), dtype=float)
    gap1=np.empty(len(data), dtype=float)
    acc1=np.empty(len(data), dtype=float)

    # export data gaps to numpy array to  limit gap to values >=GAP_VAL
    # !!! by reference, also original data['gap[m]'] affected by mainpul gapData
    
    gapData=data['gap[m]'].to_numpy()
    for i in range(0,len(data)):
        gapData[i]=max(GAP_MIN,gapData[i])
                       
    #print(f'FCdata=',FCdata)
    #print(f'gapData=',gapData)

    # initialisation
    
    v1[0]=data.loc[0,'vx[m/s]'] 
    gap1[0]=gapData[0] #max(GAP_MIN,data.loc[0,'gap[m]'])
    acc1[0]=CF.calcAccLong(gap1[0],v1[0],data.loc[0,'lead_vx']) 

    # simulation
    
    for i in range(1,len(data)):
        #if i<1181:
        if False:
            print(f'\nsim: time step i=',i,' gap1[i-1]=',gap1[i-1],' v1[i-1]=',v1[i-1])
            
        v1[i]=v1[i-1]+acc1[i-1]*dt 
        v_lead=0.5*(data.loc[i-1,'lead_vx']+data.loc[i,'lead_vx'])
        gap1[i]=gap1[i-1]+(v_lead-0.5*(v1[i]+v1[i-1]))*dt

        # if estimated speed negative, assume a stop and no further decel
        # (before possible gap reset because gap reset is dominating the actions)
        
        if v1[i]<-1e-6:  # then acc1 strictly<0
            v1[i]=0
            gap1[i]=gap1[i-1]+v_lead*dt -(-0.5*v1[i-1]**2/acc1[i-1])

        # reset gap if new leader
        
        if data.loc[i,'leadID']!=data.loc[i-1,'leadID']:
            #v1[i]=data.loc[i,'vx[m/s]']  #!! No speed reset!
            gap1[i]=gapData[i] #data.loc[i,'gap[m]']
            #print(f' new leader, i={i}: reset gap1[i]={gap1[i]}, unchanged v={v1[i]}')



        acc1[i]=CF.calcAccLong(gap1[i],v1[i],data.loc[i,'lead_vx'])

        if False:
        #if i<10:
            print(f'i={i} sLast={gap1[i-1]} vLast={v1[i-1]} vlLast={data.loc[i-1,"lead_vx"]} acc1[i-1]={acc1[i-1]}')

    # re-convert arrays to dataframe to be consistent
    # (this one-shot conversion is fast)
    
    data['v1']=v1.tolist()
    data['gap1']=gap1.tolist()
    data['acc1']=acc1.tolist()

    sse= sum((data['gap1']-data['gap[m]'])**2)
    sse1= sum((gapData-gap1)**2)
    #avg_error=sum(abs(data['gap1']-data['gap[m]']))/len(data) Ankit: measure MAD
    avg_error=np.sqrt(sse/len(data))
    #print(f'simulated CF pair with v0={v0},T={T},s0={s0},a={a},b={b}, SSE={sse}, sse1={sse1}')
    print(f'simulated CF pair with v0={v0},T={T},s0={s0},a={a},b={b}, SSE={sse}')
    return  [sse,avg_error,data]    
 



#function to load trajectory data in csv format from FCdata files
"""
@param i: subject ID

return data
"""
def singleleader(i,j):
    #path=fr"/home/mtreiber/Downloads/d8_0900_0930_road{i}_veh{j}.FCdata"  #OK
    path=fr"d8_0900_0930_road{i}_veh{j}.FCdata"
    _data=[]
    with open(path,'r') as f:

        for index,line in enumerate(f):

            if index ==8:
                as_list=line.split()
                as_list[0]='time[s]'
                header=as_list
            if index>8:
                _data.append(line.split())

    FCdata=pd.DataFrame(_data,columns=header)
    FCdata=FCdata.apply(pd.to_numeric,errors='coerce') # to convert data from object to numeric with NA
    return FCdata


###################################################################
# main
###################################################################


#load data: singleleader(roadID,vehID)

FCdata145=singleleader(2,145)
FCdata157=singleleader(2,157)
FCdata1085=singleleader(4,1085)
FCdata1377=singleleader(4,1377)

#sim

sim145=sim(v0=13.785,T=1.522,s0=0.59,a=1.670,b=5.252,FCdata=FCdata145)
#quit()
sim157=sim(v0=6.882,T=1.404,s0=0.421,a=0.969,b=4.701,FCdata=FCdata157)
sim1085=sim(v0=21.397, T=0.862, s0=1.388, a=1.445, b=3.608,FCdata=FCdata1085)
sim1377=sim(v0=9.047,T=2.160,s0=1.346,a=1.571,b=1.920,FCdata=FCdata1377)

quit()

# test 20 calculations (scan over v0)

ngrid=20
scan145_v0=np.empty(len(FCdata145), dtype=object)
for iv0 in range(0, ngrid):
    v0scan=5+iv0*1
    scan145_v0[iv0]=sim(v0=v0scan,T=1.522,s0=0.59,a=1.670,b=5.252,FCdata=FCdata145)



quit()
######################################################################


def speed_gap_plt(data):
    fig, axes=plt.subplots(2,1)
    ax=axes[0]
    ax.plot(data['time[s]'],data['vx[m/s]'],label='Data')
    ax.plot(data['time[s]'],data['v1'],label='IDM')
    ax.plot(data['time[s]'],data['lead_vx'],label='Leader')
    ax.set_ylabel('Speed')
    ax.legend()
    ax=axes[1]
    ax.plot(data['time[s]'],data['gap[m]'],label='Data')
    ax.plot(data['time[s]'],data['gap1'],label='IDM')
    ax.legend()
    ax.set_ylabel('Gap')
    ax.set_xlabel('Time')
    plt.show()

#speed_gap_plt(sim145[2])



sim_b0_b1_1=sim(v0=3.883500,T=1.280662,s0=0.59,a=1.670,b=5.252,FCdata=FCdata145)

print(f'1. sse is {sim_b0_b1_1[0]}')

sim_b0_b1_2=sim(v0=9.095485,T=1.586720,s0=0.59,a=1.670,b=5.252,FCdata=FCdata145)

print(f'2. sse is {sim_b0_b1_2[0]}')

sim_b0_b1_3=sim(v0=22.125448,T=1.618937,s0=0.59,a=1.670,b=5.252,FCdata=FCdata145)

print(f'3. sse is {sim_b0_b1_3[0]}')


sim_b2_b4_1=sim(v0=13.7852,T=1.52229,s0=0.254847,a=1.670,b=0.050000,FCdata=FCdata145)

print(f'1. sse is {sim_b2_b4_1[0]}')

sim_b2_b4_2=sim(v0=13.7852,T=1.52229,s0=0.366526,a=1.670,b=8.673333,FCdata=FCdata145)

print(f'2. sse is {sim_b2_b4_2[0]}')

sim_b2_b4_3=sim(v0=13.7852,T=1.52229,s0=0.411198,a=1.670,b=5.025000,FCdata=FCdata145)

print(f'3. sse is {sim_b2_b4_3[0]}')
