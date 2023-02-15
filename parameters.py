Written by: Dr Rea Laila Antoniou Kourounioti

import math
import numpy as np
import csv

# parameters for 
s1=			0.016 #ps1
s2=			0.0111
s3=			0.75
sel=		0.1
r1=			0.05 #pr1
r2=			0 #pr2
k=			0.18
n0i=		0.2 #proportion of silenced cells initially
Tq1=		-1
Tq2=		18
nT1=		11.5
nT2=		15
kFLC=       0.04

# constant parameters
dn=			32
g1=         0.4
growthe=	0.22
f=  		2.77 #pf
l=          2.77
# VIN3 parameters
vV=4
warm=15
dV=18
dv=0
A1=0.75
BT=17
dB=0.009
CT1=8
CT2=15.4
C1=0.0315
C2=0.03
D1=2.05
sv=vV*dV

# related to light and time of day
lightON=10
lightOFF=18
moonSt=16
moonEnd=18
sunrise=np.array([[0,lightON,lightOFF]])
sunrise=sunrise.astype(float)

# initial conditions
tEnd=0
wi=1;# vi=[0, 0];
COLFRInv=(1-n0i)*r1/(s1+r1)
n1i=(1-n0i)*s1/(s1+r1)
init=[];#to be defined later

# example parameters and temperature conditions
tiTi=np.array([[-5,-0.01,0,28,28.01,100],[22,22,5,5,22,22]])
param=[s1,s2,s3,sel,r1,r2,k,n0i,Tq1,Tq2,nT1,nT2,kFLC, f]

print("%%%%%%%%%%%%%%%%%parameters imported successfully%%%%%%%%%%%%%%")
