# -*- coding: utf-8 -*-
import numpy as np
#plt.rcParams['text.usetex'] = False


#A STATISTICAL RE-ANALYSIS OF RESISTANCE AND PROPULSION DATA - J. Holtrop
#

#Caratteristiche acqua - da inserire
rho = 1025 #densità kg/m3
ni = 1.19E-06 #viscosità
g = 9.81

#Caratteristiche nave - da inserire
L = 112.45 #between perpendiculars[m]
B = 18.5 #[m]
T = 5.400 #m. draught [m]
S = 2685.8 #wetted surface [m^2]
A_M = 98.935 #midship section area [m^2]
W_A = 1787.6 #waterplane area
LCB = 57.423 #longitudinal center of buoyancy [m]
deltaCF = 0.0002 
nabla = 8491.03 #displacement [m^3]
S_app = 3.37 #appendages wetted area [m^2]
A_T = 0.387 #transverse immersed transom area at rest [m^2]
A_BT = 15.819 #transverse area of the bulbous bow [m^2]
h_B = 2.773 #vertical position of the centre of A_BT above keel plane (<0.6T_F) [m]
T_F = 5.400 #forward draught [m]
V_kn = 16 #ship speed [kn]
V=V_kn/1.944 #ship speed [m/s]

#Calcoli
C_B = nabla/(L*B*T) #block coeff.
C_P = nabla/(L*A_M) #prismatic coeff.
lcb = (LCB-L/2)*100/L #lon. center o buoy. forw. of 0.5L as perc. of L
C_stern = 5 #Afterbody form:
#Pram with gondola = -25
#V-shaped sections = -10
#Normal = 0
#U-shaped Hogner stern = +10
C_WP = W_A/(L*B) #waterplane area coefficient
C_M = A_M/(B*T)
Rn = V*L/ni #Reynolds
Fn=V/np.sqrt(L*g) #Froude (0.219)

C_F = 0.075/((np.log10(Rn)-2)**2) #flat plate friction ITTC57
R_F = 0.5*rho*(V**2)*S*C_F #frictional resistance ITTC57 [kN]

#FORM FACTOR
L_R = L*(1-C_P+0.06*C_P*lcb/(4*C_P-1))
c_14 = 1+0.011*C_stern
k1 = 0.93+0.487118*c_14*(B/L)**1.06806*(T/L)**0.46106*(L/L_R)**0.121563*(L**3/nabla)**0.36486*(1-C_P)**-0.604247
R_F=R_F*k1

#APPENDAGE RESISTANCE
k2 = 2.8 #appenages coeff.
R_app = 0.5*rho*V**2*S_app*k2*C_F #resistance of appendages [N]

#WAVE RESISTANCE Fn>0.5
#da fare, in futuro, magari

#WAVE RESISTANCE 0.4<Fn<0.5
#da fare, in futuro, magari

#WAVE RESISTANCE Fn<0.4
if B/L<=0.11:
    c_7 = 0.229577*(B/L)**0.33333
elif 0.11<B/L<0.25:
    c_7=B/L
else:
    c_7=0.5-0.0625*L/B
i_E=1+89*np.exp(-(L/B)**0.80856*(1-C_WP)**0.30484*(1-C_P-0.0225*lcb)**0.6367*(L_R/B)**0.34574*(100*nabla/L**3)**0.16302)
# ^ half angle of entrance [degrees]
c_1=2223105*c_7**3.78613*(T/B)**1.07961*(90-i_E)**-1.37565
c_3 = 0.56*A_BT**1.5/(B*T*(0.31*np.sqrt(A_BT)+T_F-h_B))
c_2 = np.exp(-1.89*np.sqrt(c_3))
c_5 = 1-0.8*A_T/(B*T*C_M)
d = -0.9
if C_P<0.8:
    c_16 = 8.07981*C_P-13.8673*C_P**2+6.984388*C_P**3
else:
    c_16 = 1.73014-0.7067*C_P
m_1 = 0.0140407*L/T-1.75254*nabla**(1/3)/L-4.79323*B/L-c_16
if L**3/nabla<=512:
    c_15=-1.69385
elif 512<L**3/nabla<1726.91:
    c_15=-1.69385+(L/nabla**(1/3)-8)/2.36
else:
    c_15=0
m_4 = c_15*0.4*np.exp(-0.034*Fn**-3.29)
if L/B <=12:
    lmbd = 1.44*C_P-0.03*L/B
else:
    lmbd = 1.446*C_P-0.36
R_W=c_1*c_2*c_5*nabla*rho*g*np.exp(m_1*Fn**d+m_4*np.cos(lmbd*Fn**-2)) #wave resistance [N]

#BULBOUS BOW
P_B = 0.56*np.sqrt(A_BT)/(T_F-1.5*h_B)
Fni = V/np.sqrt(g*(T_F-h_B-0.25*np.sqrt(A_BT))+0.15*V**2)
R_B = 0.11*np.exp(-3*P_B**-2)*Fni**3*A_BT**1.5*rho*g/(1+Fni**2) #bulbous bow resistance [N]

#TRANSOM
FnT = V/np.sqrt(2*g*A_T/(B+B*C_WP))
if FnT<5:
    c_6 = 0.2*(1-0.2*FnT)
else:
    c_6=0
R_TR = 0.5*rho*V**2*A_T*c_6 #transom resistance [N]

#MODEL-SHIP CORRELATION
if T_F/L <= 0.04:
    c_4=T_F/L
else:
    c_4=0.04
C_A = 0.006*(L+100)**-0.16-0.00205+0.003*np.sqrt(L/7.5)*C_B**4*c_2*(0.04-c_4)
R_A = 0.5*rho*V**2*S*C_A #model-ship correlation resistance [N]

#TOTAL RESISTANCE
R_TOT = R_F+R_app+R_W+R_B+R_TR+R_A #[N]

#EFFECTIVE POWER
P_E=R_TOT*V #[W]

#RISULTATI
print(V_kn)
print(Fn)
print(R_F)
print(R_W)
print(R_B)
print(R_TR)
print(R_A)
print(R_TOT)
print(P_E)