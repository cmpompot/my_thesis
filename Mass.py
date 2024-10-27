
import scipy.constants
import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
#prwto plot

R1=[1.044*10**4,2.62*10**4,6.59*10**4,1.65*10**5,4.17*10**5,1.04*10**6,
2.62*10**6,6.58*10**6,8.29*10**6,1.655*10**7,3.302*10**7,6.589*10**7,
1.315*10**8,2.624*10**8,3.304*10**8,5.237*10**8,8.301*10**8,1.045*10**9,
1.316*10**9,1.657*10**9,2.626*10**9,4.164*10**9,6.601*10**9,8.312*10**9,
1.046*10**10,1.318*10**10,1.659*10**10,2.090*10**10,2.631*10**10,3.313*10**10,
4.172*10**10,5.254*10**10,6.617*10**10,8.332*10**10,1.049*10**11,1.322*10**11,
1.664*10**11,1.844*10**11,2.096*10**11,2.640*10**11,3.325*10**11,4.188*10**11,4.299*10**11]

P1=[9.744*10**18,4.968*10**19,2.431*10**20,1.151*10**21,5.266*10**21,2.318*10**22,
9.755*10**22,3.911*10**23,5.259*10**23,1.435*10**24,3.833*10**24,1.006*10**25,
2.604*10**25,6.676*10**25,8.738*10**25,1.629*10**26,3.029*10**26,4.129*10**26,
5.036*10**26,6.860*10**26,1.272*10**27,2.356*10**27,4.362*10**27,5.662*10**27,
7.702*10**27,1.048*10**28,1.425*10**28,1.938*10**28,2.503*10**28,3.404*10**28,
4.628*10**28,5.949*10**28,8.089*10**28,1.100*10**29,1.495*10**29,2.033*10**29,
2.597*10**29,2.892*10**29,3.290*10**29,4.473*10**29,5.816*10**29,7.538*10**29,7.805*10**29]
P11=[p*10**(-4) for p in P1]
z=[p/(r*scipy.constants.c**2) for r,p in zip(R1,P11)]
t=np.linspace(0.001,0.6,1000)
z0=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-(math.sin(x)*math.cos(x))
	z0.append((a/(b-c/d))-1)
x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
#M1=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R1,x)]
M1=[7.57292*10**(8)*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R1,x)]	
 #7.57292*10^8*(1/(ρ))^(0.5) * (sin(χ))^3

zc=max(z)
u=[scipy.constants.c,0.75*scipy.constants.c,0.57736*scipy.constants.c,0.5*scipy.constants.c,0.25*scipy.constants.c]
gammac=[((zc+1)/zc)*(speed/scipy.constants.c)**2 for speed in u]
#Mc=[5.027*10**(-31)*0.5*(3/(8*3.14))**0.5*(scipy.constants.c**2/(6.67*10**(-11)))**1.5*((1-(zc+1)/gammac)/(max(R1)/10**(-3)-max(P1)*10**(-1)/U**2))**0.5*math.sin(x[z.index(zc)])**3 for U in u]


#0.05316490009919698
Mc=[]
for U,gam in zip(u,gammac):
	q=(3/(8*3.14))**0.5
	w=(scipy.constants.c**2/(6.67*10**(-11)))**1.5
	e=(1-(1+1)/gam)
	r=(5*10**14/10**(-3)-7*10**14*33**(-1)/U**2)
	Mc.append(5.027*10**(-31)*0.5*q*w*(e/r)**0.5*math.sin(x[z.index(zc)]))



plt.plot(x[M1.index(max(M1))],max(M1),"bo",label="Current Calculation")
plt.plot(x,M1,"b-")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="upper right")
plt.grid("True")




#2o plot
R2=[7.86,7.90,8.15,11.6,16.4,45.1,212,1150]
P2=[0.9*10**9,1.01*10**10,1.01*10**11,1.21*10**12,1.40*10**13,1.70*10**14,5.82*10**15,1.90*10**17]
P22=[p*10**(-4) for p in P2]
z=[p/(r*(scipy.constants.c**2)) for r,p in zip(R2,P22)]
t=np.linspace(0.0002,0.14,1000000)
z0=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-(math.sin(x)*math.cos(x))
	z0.append((a/(b-c/d))-1)
x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
M2=[75.7292 * (10**(14) /r)**(1/2) * (math.sin(w))**3 for r,w in zip(R2,z)]	
plt.title(" Feynman-Metropolis-Teller")
plt.grid("True")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.plot(x[M2.index(max(M2))],max(M2),"ro",label="Maximum mass at 1.3*10$^{-13}$ solar masses")
plt.plot(x,M2,"r--",)
plt.legend(loc="upper left")
plt.show()
"""



"""
#3o plot
t=np.linspace(0.01,1.4,10000)
z0=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-(math.sin(x)*math.cos(x))
	z0.append((a/(b-c/d))-1)
R3=[4.460*10**11,5.228*10**11,6.610*10**11,7.964*10**11,9.728*10**11,
1.196*10**12,1.471*10**12,1.805*10**12,2.202*10**12,2.930*10**12,3.833*10**12,
4.933*10**12,6.248*10**12,7.801*10**12,9.611*10**12,1.246*10**13,1.496*10**13,
1.778*10**13,2.210*10**13,2.988*10**13,3.767*10**13,5.081*10**13,6.193*10**13,
7.732*10**13,9.826*10**13,1.262*10**14,1.586*10**14,2.004*10**14,2.520*10**14,
2.761*10**14,3.085*10**14,3.433*10**14,3.885*10**14,4.636*10**14,5.094*10**14]

P3=[7.890*10**29,8.352*10**29,9.098*10**29,9.831*10**29,1.083*10**30,1.218*10**30,
1.399*10**30,1.638*10**30,1.950*10**30,2.592*10**30,3.506*10**30,4.771*10**30,
6.481*10**30,8.748*10**30,1.170*10**31,1.695*10**31,2.209*10**31,2.848*10**31,
3.931*10**31,6.178*10**31,8.774*10**31,1.386*10**32,1.882*10**32,2.662*10**32,
3.897*10**32,5.861*10**32,8.595*10**32,1.286*10**33,1.900*10**33,2.242*10**33,
2.751*10**33,3.369*10**33,4.286*10**33,6.103*10**33,7.391*10**33]
P33=[p*10**(-4) for p in P3]
z=[p/(r*(scipy.constants.c**2)) for r,p in zip(R3,P33)]
x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
M3=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R3,x)]
#M3=[75.7292 * (10**(14) /r)**(1/2) * (math.sin(w))**3 for r,w in zip(R3,x)]	
#plt.title("Baym-Bethe-Pethick")
plt.grid("True")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.plot(x[M3.index(max(M3))],max(M3),"go",label="Baym-Bethe-Pethick")
plt.plot(x,M3,"g:")
plt.legend(loc="best")
plt.legend()

print(z)




t=np.linspace(0.1,1.4,1000)
z0=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-(math.sin(x)*math.cos(x))
	z0.append((a/(b-c/d))-1)
R4=[6.968*10**14,8.797*10**14,1.066*10**15,1.448*10**15,1.854*10**15,2.399*10**15,3.125*10**15,4.086*10**15,5.369*10**15,7.740*10**15,1.129*10**16,1.729*10**16,2.607*10**16,3.728*10**16,5.112*10**16]
P4=[2.72*10**34,3.789*10**34,5.478*10**34,1.11*10**35,2.153*10**35,4.108*10**35,7.522*10**35,1.36*10**36,2.379*10**36,4.616*10**36,8.437*10**36,1.525*10**37,2.574*10**37,3.943*10**37,5.645*10**37]
P44=[p*10**(-4) for p in P4]

z=[p/(r*scipy.constants.c**2) for r,p in zip(R4,P44)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
x=[]
z=[5*10**34*10**(-4)/(1*10**15*9*10**16)]
print(z)
x.append(np.interp(z,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
print(x)
input()

#M4=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R4,x)]
M4=[753460178.06721*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R4,x)]	
print(x[M4.index(max(M4))])
print(R4[M4.index(max(M4))])
print(max(M4))
print(753460178.06721*(1/(3125000000000000))**(0.5) * (math.sin(1.05))**3)
#M44=[75.7292 * (10**(14) /(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R4,x)]
#plt.title("Panharipande for hyperionic matter")
plt.plot(x[M4.index(max(M4))],max(M4),"mo",label="Panharipande for hyperionic matter")
plt.plot(x,M4,"m-")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="best")
plt.grid("True")
plt.show()




R5=[7.004*10**14,8.877*10**14,1.082*10**15,1.498*10**15,1.955*10**15,2.597*10**15,3.483*10**15,4.695*10**15,6.342*10**15,9.392*10**15,1.391*10**16,2.135*10**16,3.197*10**16,4.518*10**16,6.115*10**16]
P5=[3.235*10**34,5.958*10**34,9.903*10**34,2.188*10**35,4.029*10**35,7.355*10**35,1.295*10**36,2.195*10**36,3.585*10**36,6.389*10**36,1.086*10**37,1.862*10**37,3.004*10**37,4.447*10**37,6.208*10**37]
P55=[p*10**(-4) for p in P5]

z=[p/(r*scipy.constants.c**2) for r,p in zip(R5,P55)]
x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
M5=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R5,x)]
#M5=[75.7292 * (10**(14) /(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R5,x)]
plt.title("Mass from different equations of state")
plt.plot(x[M5.index(max(M5))],max(M5),"co",label="Panharipande for Pure-Neutron Matter")
plt.plot(x,M5,"c-")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="best")
plt.grid("True")
plt.show()



t=np.linspace(0.01,1.4,10000)
z0=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-(math.sin(x)*math.cos(x))
	z0.append((a/(b-c/d))-1)
	
P21=([2.163,3.687,5.822,8.681,12.38,17.04,22.79,29.74,38.02,47.76,59.09,72.14,
87.05,103.9,123,144.3,167.9,194.2,223.1,254.8,289.5,327.4,368.5,413.0,455.3,
471.1,486.9,502.8,518.6,534.5,550.3,566.2,582.1])

E21=[167.8,198.8,230.1,261.8,293.8,326.2,359.2,392.7,426.8,461.6,497.1,533.3,
570.5,608.5,647.5,687.6,728.7,771.0,814.5,859.3,905.5,953.1,1002.0,1053.0,
1105,1158,1212,1266,1321,1376,1431,1487,1544]

z=[p/e for e,p in zip(E21,P21)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R21=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E21]
M21=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R21,x)]
#M=[75.7292 * (10**(14) /(r))**(1/2) * (math.sin(w))**3 for r,w in zip(R,x)]	
plt.plot(x[M21.index(max(M21))],max(M21),"ro",label="The representative equation of state 1")
plt.plot(x,M21,"r--")
plt.title(" Calculation from Baym,Pethit and Sutherland")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="upper left")
plt.grid("True")


P22=[3.542,6.934,12.33,20.39,31.89,47.7,68.76,96.16,129.7,141.2,153,164.8,176.8,189.0,201.2,213.6,226.1,238.7]
E22=[168.5,200,232.3,265.3,299.6,335.2,372.6,412.1,454.1,497.7,542.2,587.4,633.4,680,727.3,775.2,823.8,872.9]
z=[p/e for e,p in zip(E22,P22)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R22=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E22]
M22=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R22,x)]
plt.plot(x[M22.index(max(M22))],max(M22),"bo",label="The representative equation of state $2$")
plt.plot(x,M22,"b-")
plt.title(" Calculation from Baym,Pethit and Sutherland")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="upper left")
plt.grid("True")

P23=[3.542,12.13,34.81,42.24,49.44,56.96,64.79,72.9,81.29,89.94,98.84,108,117.4,127,136.8]
E23=[168.5,200.4,234.5,270.9,308.1,346.2,384.9,424.4,464.5,505.2,546.5,588.5,631,674,717.5]
z=[p/e for e,p in zip(E23,P23)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R23=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E23]
M23=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R23,x)]
plt.plot(x[M23.index(max(M23))],max(M23),"mo",label="The representative equation of state $3$")
plt.plot(x,M23,"m-")
plt.title(" Kurkela Calculations")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.legend(loc="upper left")
plt.grid("True")
plt.show()



P31=[0.447,0.7162,0.9094,1.154,1.464,1.851,2.163,2.465,2.78,3.106,3.445,3.795,4.157,4.529,4.911,5.304,5.707,6.119,6.541,6.972,7.413,9.379,11.76,14.63,18.06,22.13,26.94,32.6,39.21]
E31=[87.9,108.2,119.5,131.5,144.4,158,167.8,183.3,198.8,214.3,229.8,245.4,261,276.6,292.2,307.9,323.5,339.2,354.9,370.6,386.4,402.2,418,434,450.1,466.3,482.7,499.2,515.9]
z=[p/e for e,p in zip(E31,P31)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R31=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E31]
M31=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R31,x)]
plt.plot(x[M31.index(max(M31))],max(M31),"bo",label="Soft")
plt.plot(x,M31,"b-")
plt.title("Hebeler, table 5")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")

plt.grid("True")



P32=[0.447,0.7162,0.9094,1.154,1.464,1.851,2.163,3.064,4.220,5.677,7.481,9.684,12.34,15.51,19.25,23.64,28.73,34.61,41.35,49.02,57.72,67.52,78.53,90.82,104.5,119.6,132,145.2,159.3]
E32=[87.9,108.2,119.5,131.5,144.4,158,167.8,183.3,198.9,214.6,230.3,246.3,262.3,278.6,295,311.7,328.6,345.7,363.2,380.9,399,417.4,436.3,455.6,475.3,495.6,516.3,537.4,559]
z=[p/e for e,p in zip(E32,P32)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R32=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E32]
M32=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R32,x)]
plt.plot(x[M32.index(max(M32))],max(M32),"ro",label="Intermediate")
plt.plot(x,M32,"r-")


P33=[0.696,1.15,1.473,1.88,2.392,3.028,3.542,3.542,5.24,7.512,10.48,14.3,20.39,28.47,38.98,52.49,69.59,80.56,92.63,105.8,120.2,135.9,152.9,171.2,190.9,212.1,234.8,259.1,285,312.6]
E33=[87.99,108.4,119.7,131.9,144.8,158.7,168.5,184.2,200.1,216.2,232.5,249.2,266.3,283.9,302.2,321.3,341.1,341.4,382.4,403.9,426,448.9,472.3,496.5,521.5,547.1,573.6,600.9,629]
z=[p/e for e,p in zip(E33,P33)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R33=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E33]
M33=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R33,x)]
plt.plot(x[M33.index(max(M33))],max(M33),"go",label="Stiff")
plt.plot(x,M33,"g:")
plt.title("Hebeler, table 5")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")

plt.grid("True")
plt.legend(loc="upper left")
plt.show()

print(z)
P34=[46.9,55.81,66.09,77.9,91.41,106.8,124.3,133.9,143.9,154.5,165.5,177,189.1,201.7,214.9,228.6,242.8,257.7,273.3,289.2,305.9,323.2,341.2,359.8,379.1,399,419.7,441.1,463.1,485.8,509.3,533.6,558.6,584.4,611,638.4,666.5]
E34=[532.8,550,567.5,585.2,603.3,621.8,640.7,659.9,679.4,699.2,719.1,739.4,759.9,780.6,801.7,823,844.6,866.5,888.7,911.2,934,957.1,980.6,1004,1028,1052,1077,1102,1128,1154,1180,1206,1233,1261,1289,1317,1345]
z=[p/e for e,p in zip(E34,P34)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R34=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E34]
M34=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R34,x)]
plt.plot(x[M34.index(max(M34))],max(M34),"bo",label="Soft")
plt.plot(x,M34,"b-")
plt.title("Hebeler, table 6")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.grid("True")
plt.legend(loc="upper left")

P35=[174.2,190,206.8,224.5,243.2,262.9,283.7,305.5,328.4,352.4,377.6,403.9,426.7,450.3,474.7,499.8,525.6,552.3,579.8,608.1,637.2]
E35=[581,603.4,626.3,649.7,673.6,698,722.9,748.3,774.3,800.8,827.9,855.6,883.8,912.6,941.8,971.6,1002,1033,1064,1096,1128]
z=[p/e for e,p in zip(E35,P35)]

x=[]
for k in z:
	x.append(np.interp(k,z0[:z0.index(max(z0))+1],t[:z0.index(max(z0))+1]))
R35=[e*10**(-9)*1.602177*10**39/scipy.constants.c**2 for e in E35]
M35=[145557930.38776*(1/(r))**(0.5) * (math.sin(w))**3 for r,w in zip(R35,x)]
plt.plot(x[M35.index(max(M35))],max(M35),"go",label="Intermediate")
plt.plot(x,M35,"g:")
plt.title("Hebeler, table 6")
plt.xlabel("χ(radians)")
plt.ylabel("M $(M_o)$")
plt.grid("True")
plt.legend(loc="upper left")
plt.legend(loc="upper left")
plt.ylim([0,1])
plt.show()
