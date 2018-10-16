# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 23:22:19 2018

@author: Konrad Urbański
"""
import numpy as np
import numpy.linalg as npl

class regP:
    def __init__(self):
        pass

    def calc(self,ref,meas,P,ogr):
        wy = (P*(ref-meas))
        if wy>ogr:
            wy=ogr
        if wy<-ogr:
            wy=-ogr
        return wy

    def version(self):
        print("v0.1")



class regPIpar:
    #pam=0
    def __init__(self):
        pass
        self.pam=0.0


    def calc(self,ref,meas,KP,KI,ogrP,ogrN,delta):
        global pam
        uchyb = (ref-meas)
        wyP = uchyb*KP
        wyI = uchyb*KI*delta + self.pam
        wy = wyP+wyI
        
        if wyI>ogrP:
            wyI=ogrP
        if wyI<ogrN:
            wyI=ogrN
        
        if wy>ogrP:
            wy=ogrP
        if wy<ogrN:
            wy=ogrN

        self.pam = wyI
        
        return wy


    def version(self):
        print("v0.11")


class timelag:
    #import numpy.linalg as npl
    def __init__(self):
        self.xold=0.0


    def calc(self,u,k,T,Ts):
        global xold
        A=-1.0/T
        B=k/T

        I=np.eye(1)
        Ad=npl.inv(I-Ts*A)
        Bd=np.dot(Ad,Ts*B)
        x=np.dot(Ad,self.xold)+np.dot(Bd,u)
    
        self.xold=x
    
        return x


#########################
class regPIDpar:
    #pam=0
    def __init__(self):
        #pass
        self.pam=0.0
        self.pamUchyb=0.0


    def calc(self,ref,meas,KP,KI,KD,ogrP,ogrN,delta):
        global pam
        uchyb = (ref-meas)
        wyP = uchyb*KP
        wyI = uchyb*KI*delta + self.pam
        wyD = ((uchyb-self.pamUchyb)/delta)*KD
        
        wy = wyP+wyI+wyD
        
        if wyI>ogrP:
            wyI=ogrP
        if wyI<ogrN:
            wyI=ogrN
        
        if wy>ogrP:
            wy=ogrP
        if wy<ogrN:
            wy=ogrN

        self.pam = wyI
        self.pamUchyb = uchyb
        
        return wy


    def version(self):
        print("v0.11")
#####################################



class gener01:
    #gen. prostokąta t1 czas1, t2 okres cyklu, w1 w2  wartośc przed i po skoku
    def __init__(self):
        pass
        #self.pam=0

    def out(self,t1,t2,w1,w2,t):
        #global pam
        
        tt=np.mod(t,t2)
        
        if tt<t1:
            wy=w1
        else:
            wy=w2

        return wy

    def version(self):
        print("v0.1")

class quality:
    def __init__(self):#,size):
        pass
        #self.reg=np.zeros((1,size))#na razie nie jest liczone
        self.quality=0.0
        self.ii=0
        self.outwindow=0.0#quality@period
        #self.outactual=0.0#na razie nie jest liczony
        self.told=-0.000001
        self.newdata=0
        
    def calc(self,inp,t,period,method):
        #input, actual time, time window width, method
        #of input preparation: 'sse', 'abs'
        global reg,ii,outwindow,outactual,told,newdata
        
        if method=='abs':
            data=np.abs(inp)
        elif method=='sse':
            data=inp*inp
        #print(data)
        tmini=np.mod(t,period)
        #print(tmini,period)
        if tmini>self.told:
            self.quality+=data
            self.newdata=0
            #print(self.quality)
        else:
            #print('teraz')
            self.newdata=1
            self.outwindow=self.quality
            #print(self.quality)
            self.quality=0.0
        
        self.told=tmini
        #outwindow quality index updated every period
        #actual value of quality index
        #synchronisation ready or not (1 or 0)
        return self.outwindow,self.quality,self.newdata




class KF:
    #assumed constant calculation perod by digital ABCD
    def __init__(self,RX,RY,QQ,RR):
        self.xhat = np.zeros((RX,1))
        self.P=np.random.rand(RX,RY)
        self.Q=np.eye(RX,RY)*QQ
        self.R=np.eye(RX,RY)*RR
        self.I=np.eye(RX,RY)

    def calc(self,Ad,Bd,Cd,u,meas):
        global xhat,P,Q,R,I

        #KALMAN START
        self.xhat=Ad@self.xhat+Bd@u
        self.P=Ad@self.P@Ad.T+self.Q
        S=Cd@self.P@Cd.T+self.R
        S2=npl.inv(S)
        K=self.P@Cd.T@S2

        resid=meas-Cd@self.xhat
        self.xhat=self.xhat+K@resid
    
        self.P=(self.I-K@Cd)@self.P
        #KALMAN END
        return self.xhat

    def version(self):
        print("v0.11")
        print("works properly for constant calculation period")

