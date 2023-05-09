# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex'] = False


#A STATISTICAL RE-ANALYSIS OF RESISTANCE AND PROPULSION DATA - J. Holtrop
#
# per navi a Fn > 0.4
#
# DA COMPLETARE




#Caratteristiche acqua
rho = 1025 #densità kg/m3
ni = 1.19E-06 #viscosità
g = 9.81

#Caratteristiche nave
L = 110
B = 18.499
T = 5.400
S = 2685.8 #sup bagnata [m2]
deltaCF = 0.0002 
C_P = 0.7345 #coeff prismatico
nabla = 8702.9 #displacement [tonnes]
lcb = 3.244 #lon center o buoy. forw. of 0.5L as perc. of L
C_stern = 0
c_14 = 1+0.011*C_stern

V_kn = 14 #velocità in nodi
V=V_kn/1.944 #vel in m/s
Rn = V*L/ni_S #Reynolds
Fn=V/np.sqrt(L*g) #Froude (0.219)

CF = 0.075/((np.log10(Rn)-2)**2) #flat plate friction ITTC57
R_F = 0.5*rho_S*(V**2)*S*CF #frictional resistance ITTC57

L_R = L(1-C_P+0.06*C_P*lcb/(4*C_P-1))
k1 = 0.93*0.487118*c_14((B/L)**1.06806)*((T/L)**0.46106)*((L/L_R)**0.121563)*(((L**3)/nabla)**0.36486)*((1-C_P)**-0.604247)

R_app = 0 #resistance of appendages

R_TR = 0 #additional pressure resistance of immersed transom stern

R_B = 0 #model-ship correlation resistance

c_17 = 6919.3*(C_M**-1.3346)*((nabla/L**3)**2.00977)*((L/B-2)**1.40692)
m_3 = -7.2035*((B/L)**0.326869)*((T/B)**0.605375)
A_BT = 0 #transverse area of the bulbous bow
h_B = #vertical position of the centre of A_BT above keel plane (<0.6T_F)
A_T = 0 #transverse immersed transom area at rest
T_F = 
c_3 = 0.56*A_BT/(B*T*(0.31*np.sqrt(A_BT)+T_F-h_B))
 = np.exp(-1.89*np.sqrt(c_3))
c_5= =  
lmbd = #IF
d = 
m_4 = 
c_15 = #if
R_WB = c_17*c_2*c_5*nabla*rho*g*np.exp(m_3*(Fn**d)+m_4*cos(lmbd*(Fn**-2)))

c_7 = #if
c_1 =
c_16 = #if
m_1 =
m_4 = 
R_WA = c_1*c_2*c_5*nabla*rho*g*np.exp(m_1*(Fn**d)+m_4*cos(lmbd*(Fn**-2))

R_WC = R_WA+(10*Fn-4)*(R_WB-R_WA)/1.5                                    
#P_E=R_TS*V_S
V_m=V/1.944 #vel in m/s

#TABELLA 1 - C_FM e C_TM in funzione di Rn
#f=open('outputDataModello.txt', 'w')
#for i in range(0,22):
#    f.write(str('%3.3f' % V_M[i]) +' & '+ str('%3.3f' % Rn_M[i])+ ' &' +str('%3.3f' % CT_M[i]) +' & '+ str('%3.3f' % CF_M[i]) +  r"\\" + '\n' )
#f.close()

#FIGURA 1
#fig=plt.figure(figsize=(8,4))
#ax=fig.add_subplot(111)
#ax.plot(Rn_M, CF_M, color='blue', label='$C_{FM}=\dfrac{0.075}{(log_{10}Rn_M-2)^2}$')
#ax.plot(Rn_M, CT_M, color='red', label='$C_{TM}$')
#ax.set_xlabel('$Rn_M$', fontsize='16')
#ax.set_ylabel('$C_{TM}$, $C_{FM}$', fontsize='16')
#plt.legend()
#plt.grid(True)
#plt.show()

#TABELLA 2 - C_FS e C_TS in funzione di Rn
#f=open('outputDataVero.txt', 'w')
#for i in range(0,22):
#    f.write(str('%3.3f' % V_S[i]) +' & '+ str('%3.3f' % Rn_S[i])+ ' &' +str('%3.4f' % CT_S[i]) +' & '+ str('%3.4f' % CF_S[i])+  r"\\" + '\n' )
#f.close()
