import math
from CoolProp.CoolProp import PropsSI

# ==========================================================
# DADOS MEDIDOS NO EXPERIMENTO
# ==========================================================
casos = [
    {
        "nome": "Sem vento",
        "V": 32.00,       # tensão [V]
        "I": 0.394,       # corrente [A]
        "Ts": 86.7,       # temperatura da barra [°C]
        "Tinf": 24.0,     # temperatura ambiente [°C]
        "U": 0.0,         # velocidade do ar [m/s]
    },
    {
        "nome": "Com vento",
        "V": 32.00,
        "I": 0.394,
        "Ts": 34.5,
        "Tinf": 24.0,
        "U": 3.00,
    },
]

# Dados fornecidos no roteiro
d = 25.4e-3             # diâmetro [m]
L = 195e-3              # comprimento [m]
epsilon = 0.60          # emissividade
g = 9.78                # gravidade [m/s²]
Patm = 924.5e2          # pressão às 17 h na estação Pampulha [Pa]
r = 1.20                # razão de pressão do soprador
R = 287.0               # constante do ar [J/(kg K)]
sigma = 5.67e-8         # constante de Stefan-Boltzmann

# Incertezas fornecidas no roteiro
uV = 0.03               # 3%
uI = 0.03               # 3%
ud = 0.1e-3             # [m]
uL = 0.5e-3             # [m]
uT = 1.0                # [K]

area = math.pi * d * L

for caso in casos:
    Ts = caso["Ts"] + 273.15
    Tinf = caso["Tinf"] + 273.15
    Tf = (Ts + Tinf) / 2
    deltaT = Ts - Tinf

    # Potência elétrica e coeficiente experimental
    q = caso["V"] * caso["I"]
    h_exp = q / (area * deltaT)

    # Coeficiente de radiação
    h_rad = epsilon * sigma * (Ts**2 + Tinf**2) * (Ts + Tinf)

    if caso["U"] == 0:
        # Convecção natural: Churchill e Chu
        P = Patm
        rho = PropsSI("Dmass", "T", Tf, "P", P, "Air")
        mu = PropsSI("VISCOSITY", "T", Tf, "P", P, "Air")
        k = PropsSI("CONDUCTIVITY", "T", Tf, "P", P, "Air")
        cp = PropsSI("Cpmass", "T", Tf, "P", P, "Air")

        nu = mu / rho
        alpha = k / (rho * cp)
        Pr = nu / alpha
        beta = 1 / Tf
        Ra = g * beta * deltaT * d**3 / (nu * alpha)

        Nu = (
            0.60
            + 0.387 * Ra ** (1 / 6)
            / (1 + (0.559 / Pr) ** (9 / 16)) ** (8 / 27)
        ) ** 2
        numero_adimensional = f"Ra = {Ra:.2e}"

    else:
        # Convecção forçada: Churchill e Bernstein
        P = r * Patm
        mu = PropsSI("VISCOSITY", "T", Tf, "P", P, "Air")
        k = PropsSI("CONDUCTIVITY", "T", Tf, "P", P, "Air")
        cp = PropsSI("Cpmass", "T", Tf, "P", P, "Air")

        rho = P / (R * Tinf)
        Pr = mu * cp / k
        Re = rho * caso["U"] * d / mu

        Nu = 0.30 + (
            0.62 * Re**0.5 * Pr ** (1 / 3)
            / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25
            * (1 + (Re / 282000) ** (5 / 8)) ** (4 / 5)
        )
        numero_adimensional = f"Re = {Re:.2e}"

    # Coeficiente teórico combinado
    h_conv = Nu * k / d
    h_teo = h_conv + h_rad

    # Propagação da incerteza de h_exp
    u_deltaT = math.sqrt(uT**2 + uT**2)
    incerteza_relativa = math.sqrt(
        uV**2
        + uI**2
        + (ud / d) ** 2
        + (uL / L) ** 2
        + (u_deltaT / deltaT) ** 2
    )
    u_h = h_exp * incerteza_relativa
    desvio = 100 * (h_exp - h_teo) / h_teo

    print(f"\n{caso['nome']}")
    print(f"Q = {q:.2f} W")
    print(numero_adimensional)
    print(f"Nu = {Nu:.2f}")
    print(f"h_exp = {h_exp:.2f} ± {u_h:.2f} W/(m² K)")
    print(f"h_conv = {h_conv:.2f} W/(m² K)")
    print(f"h_rad = {h_rad:.2f} W/(m² K)")
    print(f"h_teo = {h_teo:.2f} W/(m² K)")
    print(f"Desvio = {desvio:.1f}%")
