#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:25:07 2018

@author: Konrad Urbański
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl
import klasa_P_PI_PID_KF as reg

j=0
t1=0
tk=0.04*2
Ts=.1e-3
w=314.0
RR=reg.regPIpar()
RR2=reg.regPIpar()
TL1=reg.timelag()

tt=np.arange(t1,tk,Ts)
ww=np.zeros(len(tt))
cc=np.zeros(len(tt))
cc2=np.zeros(len(tt))
cc3=np.zeros(len(tt))


for t in tt:
    j+=1
    wz=10*np.sin(t*w)
    c=RR.calc(wz,0.0,0.0,150.0,5.0,-2.0,Ts)
    c2=RR2.calc(c,0,0,150,80,-80,Ts)
    c3=TL1.calc(c,1.0,0.001,Ts)
    ww[j-1]=wz
    cc[j-1]=c
    cc2[j-1]=c2
    cc3[j-1]=c3
    

plt.clf()
plt.plot(tt,ww,':',tt,cc,'--',tt,cc2,'-',tt,cc3,'-.')
plt.grid()
plt.show()