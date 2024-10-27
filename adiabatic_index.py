from decimal import Decimal as D
import math
import numpy as np 
import matplotlib.pyplot as plt
t=(0.01,1.2,0.01))

z=[]
for x in t:
	a=3*math.cos(x)
	b=(4.5)*math.cos(x)
	c=(math.sin(x))**3
	d=x-math.sin(x)*math.cos(x)
	z.append((a/(b-c/d))-1)
g=[]	
for zet in z:
	x=t[z.index(zet)]
	g.append((1+zet)*(1+((3*zet+1)/2)*((zet+1)*math.tan(x)**2/(6*zet)-1)))
plt.style.use(['seaborn-darkgrid'])
plt.title("Ο αδιαβατικός δείκτης $γ_c$")
plt.plot(z,g,"r-")
plt.xlabel("$P/ρc^2$")
plt.ylabel("$γ_c$")
plt.axis([0,0.35,0,4])
plt.show()





