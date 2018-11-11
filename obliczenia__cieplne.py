# coding=UTF-8
from __future__ import print_function
from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
from CoolProp.HumidAirProp import HAPropsSI
from math import sin
from matplotlib import pyplot as plt
import numpy as np

#parametry obiegu gazowego
t_1 = 10
t_4 = 1200
t_6 = 100 #temperatura spalin za parowaczem obiegu Clausiussa Rankine'a
p_1 = 101000
#Wartość sprężu regulowana
print("Wprowadź wartość sprężu: ")
k = float(input())
p_2 = k * p_1
print("Wartość ciśnienia za sprężarką: " + str(p_2) + " Pa")
eta_iGT = 0.85
eta_iGC = eta_iGT #zakładam taką samą sprawność wewnętrzą turbiny gazowej i sprężarki
eta_mGT = 0.98
eta_eGen = 0.97
eta_comb = 0.9
P_GT = 9000000
#parametry obiegu parowego
t_12 = 400
p_12 = 4000000
p_cond = 6000 #ciśnienie w skraplaczu
eta_iST = 0.75
eta_mST = 0.98
P_ST = 3500000
#parametry obiegu chłodzenia
t_21 = 20
t_22 = 22
p_20 = 500000
p_21 = 300000
p_22 = 100000
eta_iPump = 0.75
eta_mePump = 0.85
# OBLICZENIA CIEPLNE OBIEGU BRAYTONA
#entalpia przed sprężąrką
h_1 = PropsSI( "H", "T", (t_1 + 273.15), "P", p_1, "air")
#entropia przed sprężarką
s_1 = PropsSI( "S", "T", (t_1 + 273.15), "P", p_1, "air")
#entalpia powietrza za sprężarką (izentropowo)
h_2s = PropsSI( "H", "S", s_1, "P", p_2, "air")
#entalpia powietrza za sprężarką (rzeczywista)
h_2 = eta_iGC*(h_2s - h_1) + h_1
#entalpia przed turbiną gazową
h_4 = PropsSI( "H", "T", (t_4 + 273.15), "P", p_2, "air")
#entropia przed turbiną gazową
s_4 = PropsSI( "S", "T", (t_4 + 273.15), "P", p_2, "air")
#entalpia za turbiną (izentropowo)
h_5s = PropsSI( "H", "S", s_4, "P", p_1, "air")
#entalpia za turbiną (rzeczywista)
h_5 = -eta_iGT*(h_4 - h_5s) + h_4
#entalpia za parowaczem
h_6 = PropsSI("H", "T", (t_6 + 273.15), "P", p_1, "air")
#strumień masy powietrza przepływającego przez turbinę (liczone ze wzoru na moc elektryczną uzyskiwaną z części gazowej)
m_gas = P_GT/((h_4 - h_5s) * eta_iGT * eta_mGT * eta_eGen - (h_2s - h_1)*eta_iGC * eta_mGT)
print("Strumień masy spalin przepływajacych przez turbinę gazową: " + str(m_gas) + " kg/s")
#Ciepło oddawane do parowacza obiegu Claussiusa - Rankine'a
Q_56 = m_gas*(h_5 - h_6)
#zakładam, że wymiana ciepła w parowaczu obiegu Rankine'a jest idealna
#Q_5_6 = Q_11_12
#m_gas*(h_5 - h_6) = m_steam*(h_12 - h_11)
#OBLICZENIA CIEPLNE OBIEGU CLAUSSIUSA - RANKIENE'A
#enatlapia przed turbiną parową
h_12 = PropsSI("H", "T", t_12 + 273.15, "P", p_12, "water")
#entropia przed turbiną parową
s_12 = PropsSI( "S", "T", (t_12 + 273.15), "P", p_12, "water")
#entalpia pary za turbiną (izentropowo)
h_13s = PropsSI( "H", "S", s_12, "P", p_cond, "water")
#entalpia pary za turbiną parową (rzeczywista)
h_13 = -eta_iST*(h_12 - h_13s) + h_12
#strumień masy pary przepływającego przez turbinę parową (liczone ze wzoru na moc turbiny)
m_steam = P_ST/((h_12 - h_13s) * eta_iST * eta_mST * eta_eGen)
print("Strumień pary przepływającej przez turbinę parową: " + str(m_steam) + " kg/s")
#entalpia za skraplaczem 100% woda, stan nasycony
h_10 = PropsSI("H", "P", p_cond, "Q", 0, "water")
#przyrost entalpii wywołąny przez pompę wody zasilającej
v_10 = (PropsSI("D", "P", p_cond, "Q", 0, "water"))**(-1)
delta_h_1011 = v_10*(p_12 - p_cond)
#entalpia za pompą wody zasilającej
h_11 = h_10 + delta_h_1011
#ciepło odbierane od pierwszego obiegu
Q_1112 = m_steam*(h_12 - h_11)
#wymagana nadwyżka ciepła do zasilenia obiegu C-R
Q_HRSG = Q_56 - Q_1112
if Q_HRSG >= 0:
    print("Wymagana nadwyżka ciepła do zasilania obiegu C-R: " + str(
        Q_HRSG) + " J.")
#OBLICZENIA CIEPLNE OBIEGU CHŁODZENIA
#entalpia wody ogrzanej
h_21 = PropsSI("H", "T", t_21 + 273.15, "P", p_21, "water")
#entalpia wody przed pompą
h_22 = PropsSI("H", "T", t_22 + 273.15, "P", p_22, "water")
#przyrost entalpii wywołany pompą wody chłodzącej
v_22 = (PropsSI("D", "T", t_22 + 273.15, "P", p_22, "water"))**(-1)
delta_h_2022 = v_22*(p_20 - p_22)
#entalpia za pompą wody chłodzącej
h_20 = h_22 + delta_h_2022
#strumień masy wody chłodzącej wyznaczam z równania
#    m_steam* (h_13 - h_10) = m_cond*(h_21 - h_20) + m_cond(h_21 - h_22)
m_cond = m_steam*(h_13 - h_10)/-(2*h_21 - h_22 - h_20)
#OBLICZENIA MOCY POMP
#pompa wody zasilającej
P_wz = delta_h_1011 * m_steam * eta_iPump * eta_mePump
#pompa wody chłodzącej
P_cw = delta_h_2022 * m_cond * eta_iPump * eta_mePump
#SPRAWNOŚĆ BLOKU GAZOWO - PAROWEGO NETTO - liczona na podstawie parametrów turbin
eta_bNet = (P_GT + P_ST) / ((h_4 - h_2)*(1/eta_comb) * m_gas)
print("Sprawność obiegu netto (przy założeniu, że ciepło oddawanie w parowaczu jest wystarczające): " + str(eta_bNet) + " %")
#SPRAWNOŚĆ BLOKU GAZOWO - PAROWEGO BRUTTO (UWZGLĘDNIA MOC POMP)
eta_bb = (P_GT + P_ST - P_wz - P_cw)/((h_4 - h_2)*(1/eta_comb)*m_gas)
print("Sprawność obiegu brutto (uwzględnia pobór energii na pompy): " + str(eta_bb) + " %")
#SPRAWNOŚĆ OBIEGU BRAYTONA
eta_B = P_GT / ((h_4 - h_2)*(1/eta_comb) * m_gas)
print("Sprawność samego obiegu braytona: " + str(eta_B) + " %")
#SPRAWNOŚĆ OBIEGU CLAUSIUSSA - RANKINE'A
eta_CR = (P_ST - P_wz - P_cw ) / ((h_4 - h_2)*(1/eta_comb) * m_gas)
print("Sprawność samego obiegu Claussiusa-Ranakine'a: " + str(eta_CR) + " %")

xs = np.arange(0, 500, 10)
ys =  np.zeros(xs.size)
for x in range(0, xs.size):
  ys[x] = (P_GT + P_ST) / ((PropsSI( "H", "T", (t_4 + 273.15), "P", x * p_1, "air") -
                      (eta_iGC * ((PropsSI( "H", "S", s_1, "P", x * p_1, "air")) - h_1) + h_1)) *
                     (P_GT/(((PropsSI( "H", "T", (t_4 + 273.15), "P", x*p_1, "air") -
                      (PropsSI( "H", "S", (PropsSI( "S", "T", (t_4 + 273.15), "P", x * p_1, "air")), "P", p_1, "air"))) * eta_iGT * eta_mGT * eta_eGen -
                     (eta_iGC * ((PropsSI( "H", "S", s_1, "P", x * p_1, "air") - h_1) * eta_iGC * eta_mGT))))))
plt.plot(xs, ys)
plt.grid(True)
plt.xlim(0, 500)
plt.xlabel(u"spręż k = p_2/p_1")
plt.ylabel(u"Sprawność obiegu Parowo - gazowego")
plt.title(u"Wykres przedstawiający wpływ sprężu na sprawność")
plt.savefig("fig1.png", dpi = 72)
plt.show()
