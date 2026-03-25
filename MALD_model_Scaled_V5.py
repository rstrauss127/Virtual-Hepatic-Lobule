## MALD model from 
##


# %We consider y(1)=A, y(2)=N, y(3)=G, y(4)=H, y(5)=Z, y(6)=S, y(7)=L, y(8)=F with: 
# % A: Total Body APAP, g
# % N: Intracellular NAPQI Concentration, mol/cell
# % G: Intracellular GSH Concentration, mol/cell
# % H: Functional Hepatocytes, cells
# % Z: Damaged Hepatocytes, cells
# % S: Serum AST concentration, IU/L
# % L: Serum ALT concentration, IU/L
# % F: Serum clotting factor concentration

import numpy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Parameters
theta = 5.0                 # L             amount of blood

Hmax = 1.6e11               # cell         max hepatocytes
day = 8.64e4                # s


deltaZ = 5.0/8.64e4         # /s            damage hepatocyte lyse rate
r = 1.0/8.64e4              # /s            hepatocyte regen rate
eta = 5.12e13/8.64e4        # cell/mol/s    hepatocyte damage rate

alpha = 6.3/8.64e4          # /s            APAP clearing rate
deltaA = 0.33/8.64e4        # /s            APAP unconjugated clearing rate
p = 0.05                    # N/A           APAP fraction oxidized

gamma = 1.6e18/8.64e4       # cell/mol/s    GSH bind to NAPQI rate
deltaG = 2.0/8.64e4         # /s            GSH decay rate
kappa = 1.375e-14/8.64e4    # mol/cell/s    GSH produced rate
deltaN = 1.0e-4/8.64e4      # /s

deltaS = 0.92/8.64e4        # /s            AST clearing rate
deltaL = 0.35/8.64e4        # /s            ALT clearing rate
betaS = 200000.0            # IU            total quantity of AST
betaL = 84800.0             # IU            total quantity of ALT
Smin = 12.0                 # IU/L          mini AST
Lmin = 9.0                  # IU/L          mini ALT
dz = 5.0/8.64e4             # N/A           Rate of ALT/AST generation

betaF = 5.0/8.64e4          # /s            clotting clearing rate
Fmin = 0.75                 # N/A           mini clothing factor


# ODE system to be solved
def mald(t,y):
    
    ## APAP
    dy0= -(alpha*y[0]) - deltaA*y[0]
    
    ##############################################
    
    ## NAPQI
    dy1= (p*alpha*y[0]) - gamma*y[1]*y[2] - deltaN*y[1]
    ## GSH
    dy2= kappa - gamma*y[1]*y[2] - deltaG*y[2]
    
    ##############################################
    
    ## Hepatocytes status
    dy3= r*y[3]*(1-(y[3]+y[4])) - eta*y[1]*y[3]
    # dy3= - eta*y[1]*y[3]
    ## Hepatocyte damage
    dy4= eta*y[1]*y[3] - deltaZ*y[4]
    ## Hepatocyte Lyse
    dy5= deltaZ*y[4] - r*y[3]*(1-(y[3]+y[4]))
    ## Hepatocyte regeneration
    dy6= r*y[3]*(1-(y[3]+y[4]))
    
    ##############################################
    
    ## AST
    dy7= (dz*y[4]*betaS/theta) - (deltaS*(y[7]-Smin))
    ## ALT
    dy8= (dz*y[4]*betaL/theta) - (deltaL*(y[8]-Lmin))
    
    ##############################################
    
    ## Clotting
    dy9= betaF*((y[3]) - y[9])
    
    return [dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9]


# Call function to run model
def RunMALD(x, tMax, v1, v2, v3):
    global kappa
    kappa = v1
    global gamma
    gamma = v2
    global r
    r = v3
    
    # solution = solve_ivp(mald, [0, tMax], x,method='LSODA')
    solution = solve_ivp(mald, [0, tMax], x,method='Radau')
    return solution


# Call function to test model
def TestMALD():
    quantity = 20/151/160000000000
    days = 4
    # Plot results
    # x = [6.7e-13, 0, 0.88311e-14, 1, 0, 0, 1.18420e+04 ,6.731e+03, 0]
    x1 = [quantity, 0, 1.0e-14, 1, 0, 0, 0, 0 ,0, 0]
    x2 = [quantity, 0, 0.8e-14, 1, 0, 0, 0, 0 ,0, 0]
    x3 = [quantity, 0, 0.6e-14, 1, 0, 0, 0, 0 ,0, 0]
    x4 = [quantity, 0, 0.5e-14, 1, 0, 0, 0, 0 ,0, 0]
    solution1 = RunMALD(x1, 8.64e4*days, 1.575e-14/8.64e4, 1.10e18/8.64e4, 1.5/8.64e4)
    solution2 = RunMALD(x2, 8.64e4*days, 1.475e-14/8.64e4, 1.3e18/8.64e4, 1.0/8.64e4)
    solution3 = RunMALD(x3, 8.64e4*days, 1.275e-14/8.64e4, 1.6e18/8.64e4, 1.0/8.64e4)
    solution4 = RunMALD(x4, 8.64e4*days, 1.175e-14/8.64e4, 1.6e18/8.64e4, 0.5/8.64e4)

    
    variables1=solution1.y
    variables2=solution2.y
    variables3=solution3.y
    variables4=solution4.y
    temps1=solution1.t
    temps2=solution2.t
    temps3=solution3.t
    temps4=solution4.t        
    
    
    fig, axs = plt.subplots(2, 2,figsize=(20, 10))
    fig.suptitle('Model of Acetaminophen-induced Liver Disease')
    v=3
    axs[0, 0].plot(temps1,variables1[v], temps2, variables2[v], temps3, variables3[v], temps4, variables4[v])
    axs[0, 0].set_title('health')
    v=5
    axs[0, 1].plot(temps1,variables1[v], temps2, variables2[v], temps3, variables3[v], temps4, variables4[v])
    axs[0, 1].set_title('Regen')
    v=7
    axs[1, 0].plot(temps1,variables1[v], temps2, variables2[v], temps3, variables3[v], temps4, variables4[v])
    axs[1, 0].set_title('Damaged')
    v=8
    axs[1, 1].plot(temps1,variables1[v], temps2, variables2[v], temps3, variables3[v], temps4, variables4[v])
    axs[1, 1].set_title('Lysed')
    for ax in axs.flat:
        ax.set(xlabel='Time (s)', ylabel='')

    fig.tight_layout()
    return

# TestMALD()
