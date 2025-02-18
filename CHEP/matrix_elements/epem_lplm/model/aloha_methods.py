from __future__ import division
from . import wavefunctions
import cmath

def FFV1_0(F1,F2,V3,COUP):
    TMP0 = (F1[2]*(F2[4]*(V3[2]+V3[5])+F2[5]*(V3[3]+1j*(V3[4])))+(F1[3]*(F2[4]*(V3[3]-1j*(V3[4]))+F2[5]*(V3[2]-V3[5]))+(F1[4]*(F2[2]*(V3[2]-V3[5])-F2[3]*(V3[3]+1j*(V3[4])))+F1[5]*(F2[2]*(-V3[3]+1j*(V3[4]))+F2[3]*(V3[2]+V3[5])))))
    vertex = COUP*-1j * TMP0
    return vertex



def FFV1P0_3(F1,F2,COUP,M3,W3):
    V3 = wavefunctions.WaveFunction(size=6)
    V3[0] = +F1[0]+F2[0]
    V3[1] = +F1[1]+F2[1]
    P3 = [-complex(V3[0]).real, -complex(V3[1]).real, -complex(V3[1]).imag, -complex(V3[0]).imag]
    denom = COUP/(P3[0]**2-P3[1]**2-P3[2]**2-P3[3]**2 - M3 * (M3 -1j* W3))
    V3[2]= denom*(-1j)*(F1[2]*F2[4]+F1[3]*F2[5]+F1[4]*F2[2]+F1[5]*F2[3])
    V3[3]= denom*(-1j)*(-F1[2]*F2[5]-F1[3]*F2[4]+F1[4]*F2[3]+F1[5]*F2[2])
    V3[4]= denom*(-1j)*(-1j*(F1[2]*F2[5]+F1[5]*F2[2])+1j*(F1[3]*F2[4]+F1[4]*F2[3]))
    V3[5]= denom*(-1j)*(-F1[2]*F2[4]-F1[5]*F2[3]+F1[3]*F2[5]+F1[4]*F2[2])
    return V3


