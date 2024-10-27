import pandas as pd
import scipy.constants
import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate

def plot(s:str,n):
	data=np.loadtxt(s,unpack = True,delimiter="	")
	R=data[0]
	P=data[1]
	z=[p/(r) for r,p in zip(R,P)]
	t=np.linspace(0.000001,1.39,100)
	#u=np.linspace(0.15,1.39,10000)
	#t = np.concatenate((t, u))
	z0=[]
	for x in t:
		a=3*math.cos(x)
		b=(4.5)*math.cos(x)
		c=(math.sin(x))**3
		d=x-(math.sin(x)*math.cos(x))
		z0.append((a/(b-c/d))-1)
	
	
	#fsdsdfdfs
	x=[]
	R=[r*1.6*10**(-10)/(scipy.constants.c**2*10**(-39)) for r in R]
	for k in z:
		x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
	#M1=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R,x)]
	#M1=[7.57292*10**(8)*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R,x)]	
	M1=[1.36*10**8*(1/r)**0.5*(math.sin(w))**3 for r,w in zip(R,x)]
	ex=x[M1.index(max(M1))]
	den=data[0][M1.index(max(M1))]*1.6*10**(-13)/(10**(-45)*scipy.constants.c**2)
	Rad=(3*scipy.constants.c**2/(8*3.14*6.67*10**(-11)*den))**0.5*math.sin(ex)*10**-3
	name=s.replace(".dat","")
	print(f"{name}:{round(Rad,2)} km")
	
	#print(f"{name}: {round(max(M1),2)} M\u2092 ")
	plt.plot(x[M1.index(max(M1))],max(M1),"bo")
	plt.plot(x[n:],M1[n:],"b-")
	plt.title(name)
	plt.xlabel("χ(radians)")
	plt.ylabel("M $(M_o)$")
	plt.grid("True")
	plt.show()
plot("BL1.dat",70)
plot("BS.dat",60)
plot("DD2-GRDF.dat",850)
plot("DH.dat",80)
plot("FSU2H.dat",60)
plot("HHJ1.dat",70)
plot("NLD.dat",80)
plot("QS576.dat",50)
plot("SkI4.dat",70)
plot("WWF1.dat",50)