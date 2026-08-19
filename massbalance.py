import numpy as np

R = 0.08206 #atm/J/k
Streams = {

}

def set_stream(Sn=str, S = 0, B = 0, A = 0, E = 0, W = 0, CO2 = 0, H2 = 0, NH4OH = 0,NH4_salts = 0):
    masses = [round(i, 2) for i in [S,B,A,E,W,CO2,H2,NH4OH,NH4_salts]]
    Streams[Sn] = masses
    return

def system_input(target_butanol):
    Temperature = 37 + 273.15
    sugar = target_butanol/0.18 #Y_btOH g/g
    water_sugar = sugar/0.06
    NH4OH = 35.046 * ( ( (sugar * 0.2 - target_butanol * 1/6) / 59.04) + ( (sugar * 0.25 - target_butanol) / 87.10) )
    water_base = 0.7 * NH4OH / 0.3
    v_broth = 72 * (((sugar + water_sugar) / 1.021) + ((NH4OH + water_base) / 0.9))
    CO2_volumetric = (sugar * 0.6 * R * Temperature) / 0.04401  # L/h
    CO2_stripping = (v_broth - CO2_volumetric) * (0.04401 / (R * Temperature) )  # kg/h
    set_stream("S1", S=sugar, W=water_sugar)
    set_stream("S2", W=water_base, NH4OH=NH4OH)
    set_stream("S19", CO2=CO2_stripping)
    R1(v_broth)
    return

def R1(v_broth):
    S1 = Streams["S1"]
    S2 = Streams["S2"]
    S19 = Streams["S19"]
    Total = sum(S1) + sum(S2) + sum(S19)
    # Totals
    BtOH_total = S1[0] * 0.3 * 0.6
    Actn_total = (BtOH_total * 10/6) * 0.3
    EtOH_total = (BtOH_total * 10/6) * 0.1
    CO2_fermentation = S1[0] * 0.6
    H2 = S1[0] * 0.0224
    CO2_dissolved = 0.033 * (v_broth/72) * 0.04401  # kg/h

    # Gas Stripping
    BAE_broth = [BtOH_total, Actn_total, EtOH_total]
    BAE_selectivities = [33.2, 12.1, 6.3]
    BAE_condensate_mass_fractions = []
    for i in range(len(BAE_broth)):
        x = BAE_broth[i]/Total
        y = (BAE_selectivities[i] * (x/(1-x)))/(1 + BAE_selectivities[i] * (x/(1-x)))
        BAE_condensate_mass_fractions.append(y)
    B = BtOH_total * 0.63
    A = B * (BAE_condensate_mass_fractions[1]/BAE_condensate_mass_fractions[0])
    E = B * (BAE_condensate_mass_fractions[2]/BAE_condensate_mass_fractions[0])
    W = (B + A + E) / sum(BAE_condensate_mass_fractions)
    set_stream("S3", B=B, A=A, E=E, W=W, H2=H2, CO2=(S19[5] + CO2_fermentation - CO2_dissolved))

    # Neutralization
    NH4_salts = ((( S1[0] * 0.2 - EtOH_total) / 59.04) * 77.083 + ((S1[0] * 0.25 - BtOH_total) / 87.10) * 105.137)
    water = S2[7] - NH4_salts

    # Waste Treatment
    BtOH_waste = BtOH_total - B
    Actn_waste = Actn_total - A
    EtOH_waste = EtOH_total - E
    sugar = S1[0] - sum([BtOH_total, Actn_total, EtOH_total, CO2_fermentation, H2])
    water = water + S1[4] + S2[4] - W
    set_stream("Waste", S=sugar, B=BtOH_waste, A=Actn_waste, E=EtOH_waste, W=water, CO2=CO2_dissolved, NH4_salts=NH4_salts)
    F1()

def F1():
    #
    S3 = Streams["S3"]
    set_stream("S4", S=S3[0], B=S3[1], A=S3[2], E=S3[3], W=S3[4], NH4OH=S3[7], NH4_salts=S3[8])
    set_stream("S18", CO2=S3[5], H2=S3[6])
    return

def NRTL(x, T, paramenters):
    return

def Dec1():
    S4 = Streams["S4"]
    Total = sum(S4)
    w_BtOH_in = S4[1]/Total
    w_Water_in = S4[4]/Total
    # Lever rule - Dec1 layer composition
    C_aq = 0.0742 # empirical - BtOH in aqueous
    C_org = 0.797 # empirical - BtOH in organic
    w_aq = (C_org - w_BtOH_in) / (C_org - C_aq)
    w_org = (w_BtOH_in - C_aq) / (C_org - C_aq)
    m_org = Total * w_org
    m_aq = Total * w_aq
    #Assumption : BtOH distribution only relies on these weight fractions, the other components use NRTL model and the salts all remain in the aqueous layer.
    #Organic Layer
    B_org = C_org * S4[1]

    #Aqueous Layer
    B_aq = C_aq * S4[1]

system_input(63.13)
Dec1()

print("kg/h   S, B, A, E, W, CO2, H2, NH4OH, Salts")
for a in Streams:
    print(f"{a}, {Streams[a]} = {sum(Streams[a])}")

print(f"mass balance: in-out {sum(Streams['S3'])  } - {sum(Streams['S4']) + sum(Streams['S18'])}")