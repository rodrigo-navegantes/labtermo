import csv
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Medições da bancada (tabela da pág. 21 do roteiro das Aulas 10 e 11).
DADOS = {
    # Corpo de cobre: temperaturas em x1, x2 e x3, medidos a partir do topo
    "T1_C": 120.0, "T2_C": 127.7, "T3_C": 133.7,
    "VAB_V": 123.3, "I_A": 2.1,
    # Circuitos de água
    "Teb_C": 96.9,                         # água em ebulição no vaso
    "Ta1_C": 27.2, "Ta2_C": 36.7,          # entrada e saída da água de resfriamento
    "Vresf_mL": 285.0, "tresf_s": 60.0,    # água de resfriamento coletada e tempo
    "Vcond_mL": 41.0, "tcond_s": 600.0,    # condensado coletado e tempo
    # Ambiente
    "Tamb_C": 27.3,
    "Patm_kPa": 91.5,                      # atmosfera padrão a ~860 m; trocar pelo valor do INMET
}

# Geometria do corpo de cobre e incertezas (tabela do roteiro).
D_MM = 30.0
X_MM = np.array([3.0, 22.0, 31.0])
U_D_MM, U_X_MM = 0.1, np.array([0.5, 0.1, 0.1])
U_T_C = 1.0                    # termopares, °C
U_REL_EL = 0.05                # VAB e I, 5 %
U_REL_VOL, U_T_S = 0.01, 5.0   # volumes, 1 %; tempos, 5 s
U_REL_P = 0.01                 # pressão atmosférica, 1 %

K_CU, U_K_CU = 395.0, 10.0     # cobre puro perto de 110 °C, W/(m K) (Incropera, Tab. A.1)
G = 9.81                                               # m/s²
CSF = {"polido": 0.0128, "riscado": 0.0068}            # água-cobre (roteiro)
N_ROH = 1.0


# ---------------------------------------------------------------- propriedades
def propriedades(Teb_C):
    """Propriedades da água saturada em Teb (SI)."""
    T = Teb_C + 273.15
    f = lambda nome, x: PropsSI(nome, "T", T, "Q", x, "Water")
    return {"rho_l": f("D", 0), "rho_v": f("D", 1), "mu_l": f("V", 0),
            "k_l": f("L", 0), "cp_l": f("C", 0), "pr_l": f("Prandtl", 0),
            "sigma": f("I", 0), "h_lv": f("H", 1) - f("H", 0)}


def p_sat(T_C):
    return PropsSI("P", "T", T_C + 273.15, "Q", 0, "Water")


def t_sat(P_kPa):
    return PropsSI("T", "P", 1e3 * P_kPa, "Q", 0, "Water") - 273.15


# ---------------------------------------------------------------- correlações
def h_rohsenow(dTe, p, csf, n=N_ROH):
    """Rohsenow, como no roteiro [W/(m² K)]."""
    return (p["mu_l"] * p["h_lv"] * np.sqrt(G * (p["rho_l"] - p["rho_v"]) / p["sigma"])
            * (p["cp_l"] * dTe ** (2 / 3) / (csf * p["h_lv"] * p["pr_l"] ** n)) ** 3)


def h_foster_zuber(dTe, dPsat, p, forma="roteiro"):
    """Foster/Zuber [W/(m² K)].

    forma="roteiro": como no slide e no EES das Aulas 10 e 11 (fator 1,22,
    k em kW/(m K), cp em kJ/(kg K), h_lv em kJ/kg, colchete ao cubo).
    forma="original": Forster e Zuber (1955), SI, fator 0,00122, sem o cubo.
    """
    if forma == "roteiro":
        B = ((p["k_l"] / 1e3) ** 0.79 * (p["cp_l"] / 1e3) ** 0.45 * p["rho_l"] ** 0.49
             / (p["sigma"] ** 0.5 * p["mu_l"] ** 0.29 * (p["h_lv"] / 1e3) ** 0.24
                * p["rho_v"] ** 0.24))
        return 1.22 * B ** 3 * dTe ** 0.24 * dPsat ** 0.75
    B = (p["k_l"] ** 0.79 * p["cp_l"] ** 0.45 * p["rho_l"] ** 0.49
         / (p["sigma"] ** 0.5 * p["mu_l"] ** 0.29 * p["h_lv"] ** 0.24 * p["rho_v"] ** 0.24))
    return 0.00122 * B * dTe ** 0.24 * dPsat ** 0.75


# ---------------------------------------------------------------- cálculo
NOMES = ["Ts (°C)", "dTe (K)", "q_el (W)", "q'' (W/m²)", "h_exp (W/m²K)",
         "h_Rohsenow_polido (W/m²K)", "h_Rohsenow_riscado (W/m²K)",
         "h_FosterZuber (W/m²K)", "h_FosterZuber_original (W/m²K)",
         "q_evap (W)", "q_evap/q_el (-)", "Tsat(Patm) (°C)", "q_resf (W)",
         "dT/dx no cobre (K/m)", "q''_cond = k_Cu dT/dx (W/m²)", "q_cond = k_Cu (dT/dx) S (W)"]


def calcular(z):
    """z = [d, x1, x2, x3, T1, T2, T3, VAB, I, Teb, Vcond, tcond, Tamb, Patm,
            Vresf, tresf, Ta1, Ta2, k_Cu] (mm, °C, V, A, mL, s, kPa, W/(m K))."""
    d, x1, x2, x3, T1, T2, T3, VAB, I, Teb, Vc, tc, Tamb, Patm, Vr, tr, Ta1, Ta2, kcu = z
    S = np.pi * (d / 1000) ** 2 / 4
    inclinacao, Ts = np.polyfit(np.r_[x1, x2, x3] / 1000, np.r_[T1, T2, T3], 1)
    dTe = Ts - Teb
    q = VAB * I
    p = propriedades(Teb)
    dPsat = p_sat(Ts) - p_sat(Teb)

    rho_c = PropsSI("D", "T", Tamb + 273.15, "P", 1e3 * Patm, "Water")
    q_evap = rho_c * Vc * 1e-6 / tc * p["h_lv"]
    Tm = (Ta1 + Ta2) / 2 + 273.15
    rho_r = PropsSI("D", "T", Tm, "P", 1e3 * Patm, "Water")
    cp_r = PropsSI("C", "T", Tm, "P", 1e3 * Patm, "Water")
    q_resf = rho_r * Vr * 1e-6 / tr * cp_r * (Ta2 - Ta1)

    return np.r_[Ts, dTe, q, q / S, q / (S * dTe),
                 h_rohsenow(dTe, p, CSF["polido"]), h_rohsenow(dTe, p, CSF["riscado"]),
                 h_foster_zuber(dTe, dPsat, p), h_foster_zuber(dTe, dPsat, p, "original"),
                 q_evap, q_evap / q, t_sat(Patm), q_resf, inclinacao, kcu * inclinacao,
                 kcu * inclinacao * S]


def incerteza_numerica(funcao, valores, incertezas, separar=False):
    """Propagação GUM por diferenças centrais; retorna U com fator de abrangência 2."""
    termos = []
    for i, ui in enumerate(incertezas):
        dz = max(abs(valores[i]) * 1e-6, ui * 1e-3, 1e-8)
        zp, zm = valores.copy(), valores.copy()
        zp[i], zm[i] = zp[i] + dz, zm[i] - dz
        termos.append(((funcao(zp) - funcao(zm)) * ui / (2 * dz)) ** 2)
    termos = np.array(termos)
    U = 2 * np.sqrt(termos.sum(axis=0))
    return (U, termos) if separar else U


# ---------------------------------------------------------------- saída LaTeX
def _br(valor, casas=1):
    """Número no padrão brasileiro, com espaço fino nos milhares."""
    s = f"{valor:,.{casas}f}".replace(",", "X").replace(".", "{,}").replace("X", "\\,")
    return s


def _par(valor, U):
    """Arredonda U a 2 algarismos significativos e o valor na mesma casa."""
    if not U > 0:
        return _br(valor, 1), _br(0, 1)
    ordem = int(np.floor(np.log10(abs(U)))) - 1
    casas = max(0, -ordem)
    passo = 10.0 ** ordem
    return _br(round(valor / passo) * passo, casas), _br(round(U / passo) * passo, casas)


def escrever_resultados_tex(z, y, U, arquivo="resultados.tex"):
    p = propriedades(DADOS["Teb_C"])
    m = {"DataPratica": "06/10/2026 e 13/10/2026"}
    # nomes de comandos LaTeX não podem ter algarismos
    for chave, nome in (("T1_C", "Tum"), ("T2_C", "Tdois"), ("T3_C", "Ttres"),
                        ("Teb_C", "Teb"), ("Ta1_C", "Taum"), ("Ta2_C", "Tadois"),
                        ("Tamb_C", "Tamb")):
        m[nome] = _br(DADOS[chave], 1)
    m.update({"VAB": _br(DADOS["VAB_V"], 1), "Corrente": _br(DADOS["I_A"], 1),
              "Vcond": _br(DADOS["Vcond_mL"], 0), "tcond": _br(DADOS["tcond_s"], 0),
              "Vresf": _br(DADOS["Vresf_mL"], 0), "tresf": _br(DADOS["tresf_s"], 0),
              "Patm": _br(DADOS["Patm_kPa"], 1)})
    nomes = ["Ts", "DTe", "Qel", "Fluxo", "Hexp", "Hpolido", "Hriscado", "Hfz", "Hfzorig",
             "Qevap", "Fracevap", "Tsatatm", "Qresf", "Gradiente", "Fluxocond", "Qcond"]
    for nome, v, u in zip(nomes, y, U):
        if nome in ("Fluxo", "Fluxocond"):
            a, b = _par(v / 1e3, u / 1e3)          # kW/m²
        elif nome == "Fracevap":
            a, b = _par(100 * v, 100 * u)          # %
        elif nome == "Gradiente":
            a, b = _par(v / 1e3, u / 1e3)          # K/mm
        else:
            a, b = _par(v, u)
        m[nome], m["U" + nome] = a, b
    for nome, idx in (("Polido", 5), ("Riscado", 6), ("Fz", 7)):
        m["Razao" + nome] = _br(y[idx] / y[4], 1)
    m["CsfAjuste"] = _br(CSF["polido"] * (y[5] / y[4]) ** (1 / 3), 3)
    m["Fracevap"] = _br(100 * y[10], 0)
    m["Fraccond"] = _br(100 * y[15] / y[2], 0)
    m["Tsatatm"] = _br(y[11], 1)
    m.update({"Rhol": _br(p["rho_l"], 1), "Rhov": _br(p["rho_v"], 3),
              "Mul": _br(1e4 * p["mu_l"], 3), "Kl": _br(p["k_l"], 3),
              "Cpl": _br(p["cp_l"], 0), "Prl": _br(p["pr_l"], 3),
              "Tensao": _br(1e3 * p["sigma"], 2), "Hlv": _br(p["h_lv"] / 1e3, 0)})
    with open(arquivo, "w", encoding="utf-8") as arq:
        arq.write("% Gerado automaticamente por codigo-05.py. Não editar à mão.\n")
        for nome, valor in m.items():
            arq.write(f"\\newcommand{{\\{nome}}}{{{valor}}}\n")
    print(f"\n{arquivo} gerado com {len(m)} macros.")


# ---------------------------------------------------------------- gráfico
def grafico(z, y, U, arquivo="ebulicao.png"):
    d, x1, x2, x3, T1, T2, T3, VAB, I, Teb = z[:10]
    Ts, dTe, fluxo = y[0], y[1], y[3]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))

    # (a) temperaturas no cobre e extrapolação até o topo (x = 0)
    br = lambda v: f"{v:.1f}".replace(".", ",")
    xx = np.linspace(0, 33, 50)
    a1.plot(xx, Ts + y[13] * xx / 1000, color="C0", lw=1.4, label="Ajuste linear")
    a1.errorbar([x1, x2, x3], [T1, T2, T3], yerr=2 * U_T_C, fmt="o", color="C0",
                capsize=3, label="Termopares $T_1$, $T_2$, $T_3$")
    a1.errorbar(0, Ts, yerr=U[0], fmt="s", mfc="white", color="C0", capsize=3,
                label=f"$T_s$ extrapolada = {br(Ts)} °C")
    a1.axhline(Teb, color="0.45", ls=":", lw=1.2, label=f"$T_{{eb}}$ = {br(Teb)} °C")
    a1.set_xlabel("Profundidade a partir do topo, $x$ (mm)")
    a1.set_ylabel("Temperatura (°C)")
    a1.set_title("(a) Temperatura no corpo de cobre", fontsize=10)
    a1.grid(alpha=.25)
    a1.legend(fontsize=8)

    # (b) curva de ebulição: correlações e ponto experimental
    p = propriedades(Teb)
    dT = np.geomspace(2, 30, 200)
    dP = np.array([p_sat(Teb + t) - p_sat(Teb) for t in dT])
    a2.loglog(dT, h_rohsenow(dT, p, CSF["polido"]) * dT, "-", color="C1",
              label="Rohsenow, cobre polido")
    a2.loglog(dT, h_rohsenow(dT, p, CSF["riscado"]) * dT, "--", color="C2",
              label="Rohsenow, cobre riscado")
    a2.loglog(dT, h_foster_zuber(dT, dP, p) * dT, "-.", color="C4",
              label="Foster/Zuber (roteiro)")
    a2.errorbar(dTe, fluxo, xerr=U[1], yerr=U[3], fmt="o", color="C0", capsize=3,
                zorder=3, label="Experimental")
    a2.set_xticks([2, 3, 5, 10, 20, 30])
    a2.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    a2.xaxis.set_minor_formatter(plt.NullFormatter())
    a2.set_xlabel(r"Excesso de temperatura, $\Delta T_e = T_s - T_{eb}$ (K)")
    a2.set_ylabel("Fluxo de calor, $q''$ (W/m²)")
    a2.set_title("(b) Curva de ebulição nucleada", fontsize=10)
    a2.grid(alpha=.25, which="both")
    a2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(arquivo, dpi=300)
    plt.show()


# ---------------------------------------------------------------- principal
def main():
    obrig = ["T1_C", "T2_C", "T3_C", "VAB_V", "I_A", "Teb_C", "Ta1_C", "Ta2_C",
             "Vresf_mL", "tresf_s", "Vcond_mL", "tcond_s", "Tamb_C", "Patm_kPa"]
    faltando = [k for k in obrig if np.isnan(DADOS[k])]
    if faltando:
        raise SystemExit("Preencha em DADOS: " + ", ".join(faltando))

    D = DADOS
    z = np.r_[D_MM, X_MM, D["T1_C"], D["T2_C"], D["T3_C"], D["VAB_V"], D["I_A"],
              D["Teb_C"], D["Vcond_mL"], D["tcond_s"], D["Tamb_C"], D["Patm_kPa"],
              D["Vresf_mL"], D["tresf_s"], D["Ta1_C"], D["Ta2_C"], K_CU]
    uz = np.r_[U_D_MM, U_X_MM, [U_T_C] * 3, U_REL_EL * D["VAB_V"], U_REL_EL * D["I_A"],
               U_T_C, U_REL_VOL * D["Vcond_mL"], U_T_S, U_T_C, U_REL_P * D["Patm_kPa"],
               U_REL_VOL * D["Vresf_mL"], U_T_S, U_T_C, U_T_C, U_K_CU]
    y = calcular(z)
    if y[1] <= 0:
        raise SystemExit("Ts extrapolada ficou abaixo de Teb: confira T1, T2, T3 e Teb.")
    U = incerteza_numerica(calcular, z, uz)

    with open("resultados_pratica5.csv", "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(["grandeza", "valor", "U (fator 2)"])
        w.writerows(zip(NOMES, np.round(y, 4), np.round(U, 4)))
    for nome, v, u in zip(NOMES, y, U):
        print(f"{nome:32s} = {v:12.5g} ± {u:.2g}")

    escrever_resultados_tex(z, y, U)
    grafico(z, y, U)


if __name__ == "__main__":
    main()
