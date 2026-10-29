"""Prática 6 — condensação interna na serpentina do condensador.

Reproduz o programa EES do roteiro das Aulas 12 e 13 e acrescenta o cálculo
pedido no relatório: o coeficiente de Dobson e Chato local, dividindo a região
de condensação em N elementos, e a média desses coeficientes.

Requer numpy, matplotlib e CoolProp (pip install CoolProp).
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Medições da bancada (tabela do roteiro; são as mesmas da Prática 5).
DADOS = {
    "VAB_V": 123.3, "I_A": 2.1,
    "Tcond_C": 96.9,                       # temperatura de ebulição/condensação
    "Ta1_C": 27.2, "Ta2_C": 36.7,          # entrada e saída da água de resfriamento
    "Vresf_mL": 285.0, "tresf_s": 60.0,    # água de resfriamento coletada e tempo
    "Vcond_mL": 41.0, "tcond_s": 600.0,    # condensado coletado e tempo
    "Tamb_C": 27.3,
}

# Geometria do condensador (slide "Avaliação da geometria do condensador").
GEO = {"d_t_i": 34.0e-3, "d_t_e": 39.5e-3,   # tubo envelope
       "d_h_i": 22.7e-3, "d_h_e": 32.1e-3,   # hélice da serpentina
       "d_s_i": 6.5e-3, "d_s_e": 8.5e-3,     # tubo da serpentina
       "L_t": 1.43}                          # comprimento da serpentina, m

K_VIDRO = 1.17     # vidro Pyrex a ~37 °C, W/(m K) (mesmo valor do EES do roteiro)
X_MEDIO = 0.5      # título médio usado no programa do roteiro
N = 20             # elementos da região de condensação (cálculo local)

U_T_C, U_REL_VOL, U_T_S, U_REL_EL = 1.0, 0.01, 5.0, 0.05


# ---------------------------------------------------------------- propriedades
def propriedades(Tcond, Tresf):
    """Água: resfriamento em Tresf; líquido e vapor saturados em Tcond (SI)."""
    w = lambda nome, T, q: PropsSI(nome, "T", T + 273.15, "Q", q, "Water")
    return {"cp_r": w("C", Tresf, 0), "k_r": w("L", Tresf, 0),
            "rho_r": w("D", Tresf, 0), "mu_r": w("V", Tresf, 0),
            "h_lv": w("H", Tcond, 1) - w("H", Tcond, 0),
            "rho_l": w("D", Tcond, 0), "rho_v": w("D", Tcond, 1),
            "mu_l": w("V", Tcond, 0), "mu_v": w("V", Tcond, 1),
            "k_l": w("L", Tcond, 0), "pr_l": w("Prandtl", Tcond, 0),
            "k_vidro": K_VIDRO}


# ---------------------------------------------------------------- correlação
def h_dobson_chato(x, m_cond, p, d):
    """Dobson e Chato, como no roteiro; x pode ser um vetor de títulos."""
    x = np.asarray(x, dtype=float)
    Re_l = 4 * m_cond * (1 - x) / (np.pi * d * p["mu_l"])
    Re_v = 4 * m_cond * x / (np.pi * d * p["mu_v"])
    chi_ll = ((1 - x) / x) ** 0.5 * (p["rho_v"] / p["rho_l"]) ** 0.5 * (p["mu_l"] / p["mu_v"]) ** 0.5
    chi_tt = ((1 - x) / x) ** 0.9 * (p["rho_v"] / p["rho_l"]) ** 0.5 * (p["mu_l"] / p["mu_v"]) ** 0.1
    turbulento = (Re_l > 2000) & (Re_v > 2000)
    if np.any((Re_l > 2000) != (Re_v > 2000)):
        print("Aviso: regime misto (líquido e vapor em regimes diferentes); usando chi_ll.")
    chi = np.where(turbulento, chi_tt, chi_ll)
    h = 0.023 * Re_l ** 0.8 * p["pr_l"] ** 0.4 * (1 + 2.22 / chi ** 0.89) * p["k_l"] / d
    return h, Re_l, Re_v, chi


# ---------------------------------------------------------------- cálculo
CHAVES = ["Q_el", "Q_resf", "Q_cond", "Q_sub", "f", "L_cond", "Re_resf", "h_resf",
          "R_tot", "R_e", "R_p", "R_i", "T_pe", "T_pi", "h_exp",
          "Re_l", "Re_v", "chi", "h_DC", "h_DC_N", "m_cond", "m_resf"]


def calcular(z, p=None, n=N):
    """z = [Tcond, Ta1, Ta2, Vresf, tresf, Vcond, tcond, VAB, I] (°C, mL, s, V, A)."""
    Tcond, Ta1, Ta2, Vr, tr, Vc, tc, VAB, I = z
    Tresf = (Ta1 + Ta2) / 2
    if p is None:
        p = propriedades(Tcond, Tresf)
    g = GEO
    r = {}
    r["m_resf"] = Vr * 1e-6 * p["rho_r"] / tr
    r["m_cond"] = Vc * 1e-6 * p["rho_r"] / tc          # como no EES do roteiro
    r["Q_resf"] = r["m_resf"] * p["cp_r"] * (Ta2 - Ta1)
    r["Q_cond"] = r["m_cond"] * p["h_lv"]
    r["Q_el"] = VAB * I
    r["Q_sub"] = r["Q_resf"] - r["Q_cond"]
    r["f"] = r["Q_cond"] / r["Q_resf"]                  # fração do calor na condensação
    r["L_cond"] = r["f"] * g["L_t"]                     # comprimento da região de condensação
    S_i = np.pi * g["d_s_i"] * r["L_cond"]
    S_e = np.pi * g["d_s_e"] * r["L_cond"]

    # h_cond experimental: resistências em série entre o vapor e a água de resfriamento
    r["Re_resf"] = 4 * r["m_resf"] / (np.pi * g["d_h_i"] * p["mu_r"])
    r["h_resf"] = 3.657 * p["k_r"] / g["d_h_i"]         # laminar, parede a temperatura constante
    r["R_e"] = 1 / (r["h_resf"] * S_e)
    r["R_p"] = np.log(g["d_s_e"] / g["d_s_i"]) / (2 * np.pi * p["k_vidro"] * r["L_cond"])
    r["R_tot"] = (Tcond - Tresf) / r["Q_cond"]
    r["R_i"] = r["R_tot"] - r["R_e"] - r["R_p"]
    r["T_pe"] = Tresf + r["Q_cond"] * r["R_e"]
    r["T_pi"] = r["T_pe"] + r["Q_cond"] * r["R_p"]
    r["h_exp"] = 1 / (r["R_i"] * S_i)

    # Dobson e Chato: título médio (como no EES) e média de N elementos locais
    h, Re_l, Re_v, chi = h_dobson_chato(X_MEDIO, r["m_cond"], p, g["d_s_i"])
    r.update(h_DC=float(h), Re_l=float(Re_l), Re_v=float(Re_v), chi=float(chi))
    r["h_DC_N"] = perfil_local(r["m_cond"], p, n)[1].mean()
    return r


def perfil_local(m_cond, p, n=N):
    """Títulos no centro de cada elemento (x cai de 1 a 0 ao longo da região) e h local."""
    x = 1 - (np.arange(n) + 0.5) / n
    return x, h_dobson_chato(x, m_cond, p, GEO["d_s_i"])[0]


def vetor(z):
    r = calcular(z)
    return np.array([r[k] for k in CHAVES])


def incerteza_numerica(funcao, valores, incertezas):
    """Propagação GUM por diferenças centrais; retorna U com fator de abrangência 2."""
    termos = []
    for i, ui in enumerate(incertezas):
        dz = max(abs(valores[i]) * 1e-6, ui * 1e-3, 1e-8)
        zp, zm = valores.copy(), valores.copy()
        zp[i], zm[i] = zp[i] + dz, zm[i] - dz
        termos.append(((funcao(zp) - funcao(zm)) * ui / (2 * dz)) ** 2)
    return 2 * np.sqrt(np.sum(termos, axis=0))


def fracao_variancia_resf():
    """Parcela (%) da variância de R_i vinda de Ta1, Ta2, Vresf e tresf."""
    z, D = z_dados(), DADOS
    uz = np.r_[U_T_C, U_T_C, U_T_C, U_REL_VOL * D["Vresf_mL"], U_T_S,
               U_REL_VOL * D["Vcond_mL"], U_T_S, U_REL_EL * D["VAB_V"], U_REL_EL * D["I_A"]]
    i, termos = CHAVES.index("R_i"), []
    for k, u in enumerate(uz):
        dz = max(abs(z[k]) * 1e-6, u * 1e-3)
        zp, zm = z.copy(), z.copy()
        zp[k], zm[k] = zp[k] + dz, zm[k] - dz
        termos.append(((vetor(zp)[i] - vetor(zm)[i]) * u / (2 * dz)) ** 2)
    t = np.array(termos)
    return 100 * t[[1, 2, 3, 4]].sum() / t.sum()


def conferir_com_ees():
    """Com as propriedades que o EES do roteiro usa, o código deve dar o mesmo resultado."""
    p_ees = {"cp_r": 4183, "k_r": 0.606, "rho_r": 995, "mu_r": 7.657e-4, "h_lv": 2.265e6,
             "rho_l": 960.6, "rho_v": 0.5382, "mu_l": 2.913e-4, "mu_v": 1.216e-5,
             "k_l": 0.664, "pr_l": 1.848, "k_vidro": 1.171}
    geo = dict(GEO)
    GEO["d_h_e"] = 32e-3
    r = calcular(z_dados(), p_ees)
    GEO.update(geo)
    print(f"Conferência com o EES: h_cond = {r['h_exp']:.1f} (EES 657,86), "
          f"h_DC = {r['h_DC']:.1f} (EES 591,9), L_cond = {r['L_cond']:.3f} (EES 1,172)")


def z_dados():
    D = DADOS
    return np.r_[D["Tcond_C"], D["Ta1_C"], D["Ta2_C"], D["Vresf_mL"], D["tresf_s"],
                 D["Vcond_mL"], D["tcond_s"], D["VAB_V"], D["I_A"]]


# ---------------------------------------------------------------- saída LaTeX
def _br(valor, casas=1):
    return f"{valor:,.{casas}f}".replace(",", "X").replace(".", "{,}").replace("X", "\\,")


def _par(valor, U):
    """Arredonda U a 2 algarismos significativos e o valor na mesma casa."""
    if not U > 0:
        return _br(valor, 1), _br(0, 1)
    ordem = int(np.floor(np.log10(abs(U)))) - 1
    passo = 10.0 ** ordem
    casas = max(0, -ordem)
    return _br(round(valor / passo) * passo, casas), _br(round(U / passo) * passo, casas)


def escrever_resultados_tex(r, U, p, arquivo="resultados.tex"):
    m = {"DataPratica": "20/10/2026 e 27/10/2026"}
    nomes = {"Q_el": "Qel", "Q_resf": "Qresf", "Q_cond": "Qcond", "Q_sub": "Qsub",
             "f": "Frac", "L_cond": "Lcond", "Re_resf": "Reresf", "h_resf": "Hresf",
             "R_tot": "Rtot", "R_e": "Rext", "R_p": "Rp", "R_i": "Ri", "T_pe": "Tpe",
             "T_pi": "Tpi", "h_exp": "Hexp", "Re_l": "Rel", "Re_v": "Rev", "chi": "Chi",
             "h_DC": "Hdc", "h_DC_N": "HdcN", "m_cond": "Mcond", "m_resf": "Mresf"}
    escala = {"m_cond": 1e3, "m_resf": 1e3}                 # g/s
    for i, k in enumerate(CHAVES):
        a, b = _par(r[k] * escala.get(k, 1), U[i] * escala.get(k, 1))
        m[nomes[k]], m["U" + nomes[k]] = a, b
    m.update({"Tcond": _br(DADOS["Tcond_C"], 1), "Taum": _br(DADOS["Ta1_C"], 1),
              "Tadois": _br(DADOS["Ta2_C"], 1), "Vresf": _br(DADOS["Vresf_mL"], 0),
              "tresf": _br(DADOS["tresf_s"], 0), "Vcond": _br(DADOS["Vcond_mL"], 0),
              "tcond": _br(DADOS["tcond_s"], 0), "VAB": _br(DADOS["VAB_V"], 1),
              "Corrente": _br(DADOS["I_A"], 1), "Tamb": _br(DADOS["Tamb_C"], 1),
              "N": f"{N}", "Kvidro": _br(K_VIDRO, 2),
              "FracEl": _br(100 * r["Q_cond"] / r["Q_el"], 0),
              "FracRe": _br(100 * r["R_e"] / r["R_tot"], 0),
              "FracRi": _br(100 * r["R_i"] / r["R_tot"], 0),
              "RazaoDC": _br(r["h_DC"] / r["h_exp"], 2),
              "RazaoDCN": _br(r["h_DC_N"] / r["h_exp"], 2),
              "DifN": _br(abs(100 * (r["h_DC_N"] / r["h_DC"] - 1)), 0),
              "UrelHexp": _br(100 * U[CHAVES.index("h_exp")] / r["h_exp"], 0),
              "Tresf": _br((DADOS["Ta1_C"] + DADOS["Ta2_C"]) / 2, 2),
              "Rhol": _br(p["rho_l"], 1), "Rhov": _br(p["rho_v"], 3),
              "Mul": _br(1e4 * p["mu_l"], 3), "Muv": _br(1e5 * p["mu_v"], 3),
              "Kl": _br(p["k_l"], 3), "Prl": _br(p["pr_l"], 3),
              "Hlv": _br(p["h_lv"] / 1e3, 0), "Kr": _br(p["k_r"], 3),
              "Cpr": _br(p["cp_r"], 0), "Rhor": _br(p["rho_r"], 1),
              "Mur": _br(1e4 * p["mu_r"], 3)})
    m.update({"HexpNom": _br(r["h_exp"], 0), "HdcNom": _br(r["h_DC"], 0),
              "HdcNNom": _br(r["h_DC_N"], 0), "FracRp": _br(100 * r["R_p"] / r["R_tot"], 0),
              "FracVarResf": _br(fracao_variancia_resf(), 0)})
    x, h = perfil_local(r["m_cond"], p)
    m.update({"HlocMin": _br(h.min(), 0), "HlocMax": _br(h.max(), 0)})
    with open(arquivo, "w", encoding="utf-8") as arq:
        arq.write("% Gerado automaticamente por codigo-06.py. Não editar à mão.\n")
        for nome, valor in m.items():
            arq.write(f"\\newcommand{{\\{nome}}}{{{valor}}}\n")
    print(f"\n{arquivo} gerado com {len(m)} macros.")


# ---------------------------------------------------------------- gráfico
def grafico(r, U, p, arquivo="condensacao.png"):
    iu = {k: i for i, k in enumerate(CHAVES)}
    br = lambda v, c=0: f"{v:,.{c}f}".replace(",", "X").replace(".", ",").replace("X", ".")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1, 1.5]})

    # (a) resistências térmicas em série
    nomes = ["Água de\nresfriamento, $R_e$", "Parede de\nvidro, $R_p$", "Condensado,\n$R_i$"]
    vals = [r["R_e"], r["R_p"], r["R_i"]]
    erros = [U[iu["R_e"]], U[iu["R_p"]], U[iu["R_i"]]]
    barras = a1.bar(nomes, vals, color=["C0", "0.6", "C1"], width=0.6,
                    yerr=erros, capsize=4, ecolor="0.25")
    for b, v in zip(barras, vals):
        a1.text(b.get_x() + b.get_width() / 2, v / 2, f"{100 * v / r['R_tot']:.0f}%".replace(".", ","),
                ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    a1.set_ylabel("Resistência térmica (K/W)")
    a1.set_title("(a) Resistências entre vapor e água de resfriamento", fontsize=10)
    a1.grid(alpha=.25, axis="y")

    # (b) coeficiente local ao longo da região de condensação
    x, h = perfil_local(r["m_cond"], p)
    zc = (1 - x) * r["L_cond"]
    a2.axhline(r["h_exp"], color="C0", lw=1.6, label=f"Experimental: {br(r['h_exp'])} W/(m²K) (valor nominal)")
    a2.plot(zc, h, "-o", color="C1", ms=3.5, lw=1.4, label=f"Dobson e Chato local ($N={N}$)")
    a2.axhline(r["h_DC_N"], color="C1", ls="--", lw=1.4,
               label=f"Média dos {N} elementos: {br(r['h_DC_N'])} W/(m²K)")
    a2.plot(r["L_cond"] / 2, r["h_DC"], "D", color="C4", ms=7, zorder=4,
            label=f"Título médio x = 0,5: {br(r['h_DC'])} W/(m²K)")
    a2.set_xlabel("Posição na região de condensação, $z$ (m)   [título $x = 1 - z/L_{cond}$]")
    a2.set_ylabel(r"$h_{cond}$ (W/(m²K))")
    a2.set_title("(b) Coeficiente ao longo da região de condensação", fontsize=10)
    a2.set_ylim(0, 1.45 * max(h.max(), r["h_exp"]))
    a2.grid(alpha=.25)
    a2.legend(fontsize=8, loc="upper left")

    virgula = plt.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    a1.yaxis.set_major_formatter(virgula)
    a2.xaxis.set_major_formatter(virgula)
    fig.tight_layout()
    fig.savefig(arquivo, dpi=300)
    plt.show()


# ---------------------------------------------------------------- principal
def main():
    conferir_com_ees()
    D = DADOS
    z = z_dados()
    uz = np.r_[U_T_C, U_T_C, U_T_C, U_REL_VOL * D["Vresf_mL"], U_T_S,
               U_REL_VOL * D["Vcond_mL"], U_T_S, U_REL_EL * D["VAB_V"], U_REL_EL * D["I_A"]]
    r = calcular(z)
    U = incerteza_numerica(vetor, z, uz)
    p = propriedades(D["Tcond_C"], (D["Ta1_C"] + D["Ta2_C"]) / 2)

    with open("resultados_pratica6.csv", "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(["grandeza", "valor", "U (fator 2)"])
        w.writerows((k, round(r[k], 6), round(u, 6)) for k, u in zip(CHAVES, U))
    for k, u in zip(CHAVES, U):
        print(f"{k:8s} = {r[k]:12.5g} ± {u:.2g}")

    escrever_resultados_tex(r, U, p)
    grafico(r, U, p)


if __name__ == "__main__":
    main()
