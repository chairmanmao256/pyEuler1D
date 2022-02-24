import pyEuler1D
import math
from matplotlib import pyplot as plt 
# Solve 1-D Euler equation using the C++ code wrapped by Cython
# This is the user-interface of the solver written by python

######################## User input ###########################

n = 1500                                   # The number of grid cells along x direction
cfl = 0.25                                 # CFL number. The time step is determined by: dt = CFL*dx
t = 0.1                                    # Time to terminate the time advancement (The equation is integrated to t)
bctype = bytes("Inf","ascii")              # Boundary condition. Currently only the infinity condition (Inf) and periodic 
                                           # condition (Periodic) are supported. Input otherwise would raise error.
FILENAME = "Output"                        # The name of the ouput file
file_num = 100                             # The number of output files

######################## Solver setup #########################

dx = 1.0/n
dtt = cfl*dx
step = math.floor(t/dtt)
writestep = math.floor(step/file_num)
pysolver = pyEuler1D.Euler1D(bctype, n)
pysolver.initialize()
pysolver.set_boundary()

######################### Solver run ##########################
for i in range(0,step):
    for j in range(0,4):
        pysolver.reconstruction_0()
        pysolver.avg_roe()
        pysolver.reconstruction_TVD()
        pysolver.avg_roe()
        pysolver.cal_roe_flux()
        dt = dtt/(4-j)
        pysolver.timeadvancement(dt,1-j)
        pysolver.set_boundary()
        #if (i+1)%writestep == 0:
        #    fnam = bytes(FILENAME+str(i)+".plt","ascii")
        #    pysolver.writefile(fnam)

fnam = bytes(FILENAME+".dat","ascii")           
pysolver.writefile(fnam)


##################### Post processing ########################
#read the outcome
filename = FILENAME +'.dat'
prof = open(filename,'r')
x = []
rho = []
u = []
p = []

lineprof = prof.readlines()
index = 0
for line in lineprof:
    coord = line.split()
    x.append(float(coord[0]))
    rho.append(float(coord[1]))
    u.append(float(coord[2]))
    p.append(float(coord[3]))


#plot the outcome
plt.subplot(1,3,1)
plt.title("Shock-tube, rho, t=" + str(t)) 
plt.xlabel("x") 
plt.ylabel("rho") 
plt.plot(x,rho) 

plt.subplot(1,3,2)
plt.title("Shock-tube, u, t=" + str(t)) 
plt.xlabel("x") 
plt.ylabel("u") 
plt.plot(x,u) 

plt.subplot(1,3,3)
plt.title("Shock-tube, p, t=" + str(t)) 
plt.xlabel("x") 
plt.ylabel("p")
plt.plot(x,p) 

#plt.show()
#plt.figure(figsize=(6, 3))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=0.7, hspace=0.2)
plt.savefig("Shock-tube_t = " + str(t)+".png", dpi=300)