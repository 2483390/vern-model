#Written by: Dr Rea Laila Antoniou Kourounioti

import numpy as np
import parameters as p
import math

def model(y,time,tiTi):

	# variables
	b=y[0]
	v=y[1]
	V=y[2]
	Hs=y[3]
	Is=y[4]
	Ns=y[5]
	Ss=y[6]
	Ps=y[7]
	#dividing cells
	H=y[8]
	I=y[9]
	N=y[10]
	S=y[11]
	P=y[12]
	#nondividing cells
	n=Hs+Is+Ns+Ss+Ps+H+I+N+S+P;
	#total cells
	r=y[13]
	#reactivation 
	FLC=y[14]
	#FLC not divided by intial condition



	# light
	if len(p.sunrise)>1:
		if time>p.lightON/24:
			if p.sunrise[math.floor(time-p.lightON/24),0]==math.floor(time-p.lightON/24):
				day=p.sunrise[math.floor(time-p.lightON/24),]
			else:
				print("sunrise position and day don't match")
				print("day "+str(p.sunrise[(p.sunrise[:,0]==math.floor(time-p.lightON/24))])+", time "+str(time))
				day=p.sunrise[(p.sunrise[:,0]==math.floor(time-p.lightON/24))][0]
		else:
			if p.sunrise[math.floor(time),0]==math.floor(time):
				day=p.sunrise[math.floor(time),]
			else:
				print("sunrise position and day don't match")
				print("day "+str(p.sunrise[(p.sunrise[:,0]==math.floor(time))])+", time "+str(time))
				day=p.sunrise[(p.sunrise[:,0]==math.floor(time))][0]
	else:
		day=p.sunrise[0]
		
	if day.size==0:
		print("problem extracting sunrise at time:")
		print(time)
		return -1
	
	p.lightON=day[1]
	p.lightOFF=day[2]
	p.moonSt=day[2]

	
	# temperature
	T=np.interp(time,tiTi[0],tiTi[1])
	nT=nightTime12h(time,p.lightON,p.lightOFF,tiTi)

	# parameters
	a=(1+(-1+p.A1)*(MAXbefore(tiTi,time,p.moonSt)>p.warm))
	c=p.C1-((T-p.CT1)/(p.CT2-p.CT1)*(T<p.CT2)*(T>p.CT1)+(T>=p.CT2))*p.C2;
	d=(p.D1+math.sin(2*math.pi*(time-(p.lightON-1)/24)))**2;
	pv=a*b*c*d;
	s1=p.s1
	s2=p.s2*(p.Tq2-T)*(T-p.Tq1)*(T>p.Tq1)*(T<p.Tq2)*V
	g=r*p.g1*math.exp(p.growthe*(T-22))
	dn=p.dn
	if (FLC<p.kFLC)&(T>p.Tq2):
		k=p.k
	else:
		k=0
	
	r1=p.r1*((nT-p.nT1)/(p.nT2-p.nT1)*(nT<p.nT2)*(nT>p.nT1)+(nT>=p.nT2))
	s3=p.s3*g
	r2=p.r2*g

	## ODEs
	out=[(T<p.BT)-p.dB*b,#dB
		pv-p.sv*v-p.dv*v,#v
		p.sv*v-p.dV*V,#V
		-(s1+s2*p.sel)*Hs+r1*Is,#Hs
		s1*Hs-(r1+s2)*Is+r2*Ps,#Is
		s2*p.sel*Hs+s2*Is-g*Ns,#Ns
		g*Ns-s3*Ss,#Ss
		s3*Ss-r2*Ps,#Ps
		dn*g*Hs-(s1+s2*p.sel)*H+r1*I,#H
		dn*r2*Ps+dn*g*Is+s1*H-(r1+s2)*I,#I
		s2*p.sel*H+s2*I,#N
		dn*g*Ns+dn*(g-s3)*Ss,#S
		dn*s3*Ss+dn*(g-r2)*Ps,#P
		-k*r,#r
		p.f*(Hs+H)/n - p.l*FLC]#FLC rate of change, deg before FLC 
		
	p.tEnd=time
	return out
	


def output(y):
	FLCexp=y[:,14]/p.COLFRInv #divides by initial value
	VIN3=y[:,2]
	return VIN3,FLCexp
#adjusted to intial condition of FLC in COLFRI


def nightTime12h(time,lightON,lightOFF,tiTi):
	midday=np.mean([lightON,lightOFF])
	morningTime=0.25+(midday-12)/24
	night=math.floor(time+1-morningTime)-1-0.5+morningTime
	if night<tiTi[0][0]:
		nT=22
	else:
		morning=math.floor(time+1-morningTime)-1+morningTime
		Tnight=np.interp(night,tiTi[0],tiTi[1])
		Tmorning=np.interp(morning,tiTi[0],tiTi[1])
		nightIn=np.searchsorted(tiTi[0],night,side='right')#where(tiTi[0]>night)[0][0]
		morningIn=np.searchsorted(tiTi[0],morning)-1#where(tiTi[0]<morning)[0][len(np.where(tiTi[0]<morning)[0])-1]
		if nightIn<morningIn:
			Tsum=np.array(tiTi[1][nightIn:(morningIn+1)]).cumsum()
			Tsum[2:]=Tsum[2:]-Tsum[:-2]
			nT=((Tnight+tiTi[1][nightIn])*(tiTi[0][nightIn]-night)/2+np.dot(Tsum[1:],np.diff(tiTi[0][nightIn:(morningIn+1)]))/2+(Tmorning+tiTi[1][morningIn])*(morning-tiTi[0][morningIn])/2)*2
		else:
			if morningIn==nightIn:
				nT=((Tnight+tiTi[1][nightIn])*(tiTi[0][nightIn]-night)/2+(Tmorning+tiTi[1][morningIn])*(morning-tiTi[0][morningIn])/2)*2
			else:
				nT=(Tmorning+Tnight)/2
	return nT


	
def MAXbefore(tiTi,t,moonEnd):
	day=math.floor(t)
	if day<tiTi[0][0]:
		maxBefore=22
	else:
		if t%1<=(moonEnd/24+10**(-10)):
			if day>0:
				thatDay=tiTi[1][(tiTi[0]>(day-1+moonEnd/24))&(tiTi[0]<=t)]
				if len(thatDay)>0:
					maxBefore=max(np.amax(np.interp([day-1+moonEnd/24,t],tiTi[0],tiTi[1])),np.amax(thatDay))
				else:
					maxBefore=np.amax(np.interp([day-1+moonEnd/24,t],tiTi[0],tiTi[1]))
			else:
				thatDay=tiTi[1][(np.array([math.floor(x) for x in tiTi[0]])==(day))&(tiTi[0]<=t)]
				if len(thatDay)>0:
					maxBefore=max(np.amax(np.interp([day,t],tiTi[0],tiTi[1])),np.amax(thatDay))
				else:
					maxBefore=np.amax(np.interp([day,t],tiTi[0],tiTi[1]))
		else:
			thatDay=tiTi[1][(tiTi[0]>(day+moonEnd/24))&(tiTi[0]<=t)]
			if len(thatDay)>0:
				maxBefore=max(np.amax(np.interp([day+moonEnd/24,t],tiTi[0],tiTi[1])),np.amax(thatDay))
			else:
				maxBefore=np.amax(np.interp([day+moonEnd/24,t],tiTi[0],tiTi[1]))
	return maxBefore
