import numpy as np

GM = 3.9859e14                  ;% meter^3/ sec^2
tdist= 4000e3                   ;% meters : Distance of target from launch point
Re=6378.4e3                     ;% meters : Radius of the earth
t_arr= ones(1571) * 1.5e5]      ;% N TODO(sparsh) update this to use numerical data depeding on t
mass0=13.1e3                    ;% kg
g=9.8                           ;% m/sec^2
r0= Re + 50e3                   ;% meter  : 65kms above he surfface of the earth
V0= 3300                        ;% m/s
gamma0=70*np.pi/180             ;% radians
phi0= (tdist - 30e3) / Re       ;% radians
rT = Re + 65e3                  ;% meters
I_sp = 250                      ;% seconds
m_dot = t_arr/(g * I_sp)        ;% Rate of change of mass

m_arr = np.array([mass0] * 1571)
time = np.arange(0,157.1,0.1)
for i in range(1,1571):
    m_arr[i] = m_arr[i] - np.trapz(m_dot[0:i],time[0:i])
