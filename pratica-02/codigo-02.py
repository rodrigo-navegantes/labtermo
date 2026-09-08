"""EMA-103 - Pratica 2: tratamento dos dados e incertezas (k=2)."""

import math as M


# Medicoes no SI: Pa, m, s, kg e graus Celsius.
x = dict(V=120.3, I=20.0, P0=0.0, T0=96.0, dH=.042, t=1800.0,
         P1=.95e5, T1=116.4, Ta=26.4, Patm=91.86e3, g=9.78,
         mv=.0558, dp=.00249, di=.247, de=.253, H=.210, eps=.60)
u = dict(V=.03*x["V"], I=.03*x["I"], P0=.05e5, T0=1.0,
         dH=.002, t=30.0, P1=.05e5, T1=1.0, Ta=1.0,
         Patm=.02*x["Patm"], g=.01, mv=.0001, dp=.0001,
         di=.001, de=.001, H=.005)


def psat(T):
    """Pressao de saturacao [Pa] pela equacao de Antoine."""
    A, B, C = ((8.07131, 1730.63, 233.426) if T <= 99 else
               (8.14019, 1810.94, 244.485))
    return 133.322368 * 10**(A - B/(C + T))


def tsat(P):
    """Inversa da equacao de Antoine [graus Celsius]."""
    A, B, C = ((8.07131, 1730.63, 233.426) if P < 1e5 else
               (8.14019, 1810.94, 244.485))
    return B/(A - M.log10(P/133.322368)) - C


def rho(T):
    """Densidade da agua liquida [kg/m3], equacao de Kell."""
    return 1000*(1-(T+288.9414)/(508929.2*(T+68.12963))*(T-3.9863)**2)


def hfg(T):
    """Calor latente de referencia [J/kg], correlacao de Watson."""
    return 2256.95e3*((647.096-T-273.15)/(647.096-373.15))**.38


def calc(z):
    DT, Ts, Ta = z["T0"]-z["Ta"], z["T0"]+273.15, z["Ta"]+273.15
    hr = z["eps"]*5.670374419e-8*(Ts+Ta)*(Ts**2+Ta**2)
    hv, hh = 1.42*(DT/z["H"])**.25, 1.32*(DT/(z["de"]/4))**.25
    Sv, Sh = M.pi*z["de"]*z["H"], M.pi*z["de"]**2/4
    Qv, Qh = (hv+hr)*Sv*DT, (hh+hr)*Sh*DT
    Qel, rhol = z["V"]*z["I"], rho(z["T0"])
    mdot = rhol*M.pi*z["di"]**2/4*z["dH"]/z["t"]
    Pvalv = z["mv"]*z["g"]/(M.pi*z["dp"]**2/4)
    return dict(Pab=z["Patm"]+z["P0"], Ptab=psat(z["T0"]),
                Pman=z["Patm"]+z["P1"], Pval=z["Patm"]+Pvalv,
                Ptf=psat(z["T1"]), Tab=tsat(z["Patm"]+z["P0"]),
                Tman=tsat(z["Patm"]+z["P1"]),
                Tval=tsat(z["Patm"]+Pvalv), Ppeso=Pvalv,
                hv=hv, hh=hh, hr=hr, Sv=Sv, Sh=Sh, Qel=Qel,
                Qv=Qv, Qh=Qh, Qper=Qv+Qh, Qev=Qel-Qv-Qh,
                rho=rhol, massa=mdot*z["t"], mdot=mdot,
                h=(Qel-Qv-Qh)/mdot, href=hfg(z["T0"]))


def U(saida):
    """Propagacao numerica; retorna a incerteza expandida 2*uc."""
    soma = 0.0
    for nome, un in u.items():
        passo, xp, xm = max(un/100, 1e-10), x.copy(), x.copy()
        xp[nome] += passo
        xm[nome] -= passo
        soma += ((calc(xp)[saida]-calc(xm)[saida])/(2*passo)*un)**2
    return 2*M.sqrt(soma)


y = calc(x)
for nome, chave, escala, unidade in [
    ("P aberta", "Pab", 1e3, "kPa"),
    ("P fechada (manometro)", "Pman", 1e3, "kPa"),
    ("P fechada (peso)", "Pval", 1e3, "kPa"),
    ("Perdas", "Qper", 1, "W"),
    ("Vazao evaporada", "mdot", 1e-3, "g/s"),
    ("Calor latente", "h", 1e6, "MJ/kg")]:
    print(f"{nome}: {y[chave]/escala:.4g} +/- {U(chave)/escala:.3g} {unidade}")
