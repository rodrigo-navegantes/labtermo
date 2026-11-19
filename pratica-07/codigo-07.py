"""Prática 7 — bancada de refrigeração por compressão de vapor (R12).

Faz os cálculos do exemplo do roteiro das Aulas 14 e 15 (programa EES com
R134a), agora com as medições da bancada:
  * pressão na entrada do evaporador pela temperatura (saturação) x manômetro;
  * estados do R12 e diagrama P-h do ciclo;
  * efeito e capacidade de refrigeração, trabalho específico e potência de
    compressão, perda de calor do compressor e COP.

Requer numpy, matplotlib e CoolProp (pip install CoolProp).
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import CoolProp.CoolProp as CP

FLUIDO = "R12"
for _f in (FLUIDO, "R134a"):
    # h = s = 0 no líquido saturado a -40 °C: mesma referência do EES e do diagrama P-h do roteiro
    CP.set_reference_state(_f, "ASHRAE")

KGF, BAR, LBM = 98.1, 100.0, 0.454    # kPa por kgf/cm², kPa por bar, kg por lbm (roteiro)

# Medições da bancada (tabelas do roteiro). Pressões manométricas; pontos do esquema do slide 3.
DADOS = {
    "T1_C": 0.0, "P1_kgf": 2.2,       # saída da válvula termostática (mesmas leituras do ponto 2)
    "T2_C": 0.0, "P2_kgf": 2.2,       # entrada do evaporador
    "T3_C": 7.1, "P3_kgf": 2.1,       # saída do evaporador
    "T4_C": 9.0, "T5_C": 14.9,        # vapor no trocador de calor: entrada e saída
    "T10_C": 37.6, "T11_C": 32.1,     # líquido no trocador de calor: entrada e saída
    "T6_C": 20.0, "P6_kgf": 1.95,     # entrada do compressor
    "T7_C": 73.2, "P7_kgf": 9.2,      # saída do compressor
    "T8_C": 62.1, "P8_kgf": 8.6,      # entrada do condensador
    "T9_C": 39.5, "P9_kgf": 8.4,      # saída do condensador (título ~ 0)
    "T12_C": 32.2, "P12_bar": 8.2,    # entrada da válvula termostática (manômetro em bar)
    "m_lbm_min": 1.3,                 # vazão de R12 no rotâmetro
    "Tamb_C": 22.1,
    "Patm_kPa": 91.5,                 # atmosfera padrão a ~860 m (BH), como nas práticas 5 e 6
    "Tsup_C": None,                   # superfície do compressor (não medida): None -> (T6 + T7)/2
}

# Compressor visto como uma caixa a x b x c (foto do roteiro). A face a x c, voltada
# para o motor elétrico, é adiabática e fica fora da área de troca.
A, B, C = 0.10, 0.17, 0.23            # m
EPS, SIGMA, G = 0.9, 5.67e-8, 9.81

# Incertezas do roteiro (tratadas como incertezas padrão; resultados com fator 2)
U_T = 1.0                 # termopares tipo K, °C
U_M = 0.1                 # rotâmetro, lbm/min
U_PB = 0.03 * 8           # manômetros de baixa: 3 % do fundo de escala de 8 kgf/cm²
U_PA = 0.03 * 35          # manômetros de alta: 3 % de 35 kgf/cm²
U_PA_BAR = 0.03 * 34      # 4º manômetro de alta, em bar: 3 % de 34 bar
U_REL_PATM = 0.01         # pressão atmosférica, 1 %


def incertezas(d):
    u = {"T2_C": U_T, "T3_C": U_T, "T6_C": U_T, "T7_C": U_T, "T12_C": U_T, "Tamb_C": U_T,
         "P2_kgf": U_PB, "P3_kgf": U_PB, "P6_kgf": U_PB, "P7_kgf": U_PA,
         "P12_bar": U_PA_BAR, "m_lbm_min": U_M, "Patm_kPa": U_REL_PATM * d["Patm_kPa"]}
    if d["Tsup_C"] is not None:
        u["Tsup_C"] = U_T
    return u


# ---------------------------------------------------------------- propriedades
def h_TP(T_C, P_kPa, fluido=FLUIDO):
    """Entalpia (kJ/kg) de líquido ou vapor (uma só fase) a partir de T e P."""
    return CP.PropsSI("H", "T", T_C + 273.15, "P", P_kPa * 1e3, fluido) / 1e3


def h_Px(P_kPa, x, fluido=FLUIDO):
    return CP.PropsSI("H", "P", P_kPa * 1e3, "Q", x, fluido) / 1e3


def T_sat(P_kPa):
    return CP.PropsSI("T", "P", P_kPa * 1e3, "Q", 0, FLUIDO) - 273.15


def P_sat(T_C):
    return CP.PropsSI("P", "T", T_C + 273.15, "Q", 0, FLUIDO) / 1e3


def t_sup(d):
    """Temperatura da carcaça do compressor: medida, se houver; senão (T6 + T7)/2."""
    return d["Tsup_C"] if d["Tsup_C"] is not None else (d["T6_C"] + d["T7_C"]) / 2


def pressoes(d):
    """Pressões absolutas (kPa): manométrica convertida + pressão atmosférica."""
    pa = d["Patm_kPa"]
    P = {k: d[f"P{k}_kgf"] * KGF + pa for k in (2, 3, 6, 7, 8, 9)}
    P[12] = d["P12_bar"] * BAR + pa
    return P


# ---------------------------------------------------------------- compressor
def perda_compressor(Tsup, Tamb, Patm_kPa):
    """Calor perdido pela carcaça (W), com h combinado de convecção natural e radiação."""
    S = 2 * (A * B + B * C + C * A) - A * C       # caixa sem a face voltada para o motor
    L = (A * B * C) ** (1 / 3)                     # aresta do cubo de mesmo volume (L³ = abc)
    dT = Tsup - Tamb
    h_conv = 1.052 * (abs(dT) / L) ** 0.25         # Holman, tabela 7-2: cubo aquecido, laminar
    Ts, Ta = Tsup + 273.15, Tamb + 273.15
    h_rad = EPS * SIGMA * (Ts**2 + Ta**2) * (Ts + Ta)
    # Rayleigh do ar no filme, para conferir a faixa laminar da correlação (1e4 a 1e9)
    Tf, p = (Ts + Ta) / 2, Patm_kPa * 1e3
    ar = lambda nome: CP.PropsSI(nome, "T", Tf, "P", p, "Air")
    nu, alfa = ar("V") / ar("D"), ar("L") / (ar("D") * ar("C"))
    Ra = G * abs(dT) * L**3 / (Tf * nu * alfa)
    return {"S": S, "L": L, "h_conv": h_conv, "h_rad": h_rad, "Q": (h_conv + h_rad) * S * dT,
            "Ra": Ra, "Tsup": Tsup}


# ---------------------------------------------------------------- ciclo
CHAVES = ["P2_man", "P2_sat", "x2", "SA3", "SR12", "h2", "h3", "h6", "h7", "m",
          "q_evap", "Q_evap", "w", "W_ad", "h_conv", "h_rad", "h_comb", "Q_comp",
          "W_comp", "COP", "COP_ad"]


def calcular(d):
    P = pressoes(d)
    r = {}
    # pressão na entrada do evaporador: manômetro (+ P_atm) e saturação na temperatura medida
    r["P2_man"] = P[2]
    r["P2_sat"] = P_sat(d["T2_C"])

    # estados do R12 (pontos 2, 3, 6, 7 e 12 do formulário)
    r["h12"] = h_TP(d["T12_C"], P[12])              # líquido subresfriado
    r["h2"] = r["h12"]                              # expansão isentálpica na válvula termostática
    hl, hv = h_Px(P[2], 0), h_Px(P[2], 1)
    r["x2"] = (r["h2"] - hl) / (hv - hl)            # título na entrada do evaporador
    r["h3"] = h_TP(d["T3_C"], P[3])                 # vapor superaquecido
    r["h6"] = h_TP(d["T6_C"], P[6])
    r["h7"] = h_TP(d["T7_C"], P[7])
    r["SA3"] = d["T3_C"] - T_sat(P[3])              # superaquecimento na saída do evaporador
    r["SR12"] = T_sat(P[12]) - d["T12_C"]           # subresfriamento na entrada da válvula

    # balanços de energia (kJ/kg e kW)
    r["m"] = d["m_lbm_min"] * LBM / 60              # kg/s
    r["q_evap"] = r["h3"] - r["h2"]                 # efeito de refrigeração
    r["Q_evap"] = r["m"] * r["q_evap"]              # capacidade de refrigeração
    r["w"] = r["h7"] - r["h6"]                      # trabalho específico de compressão
    r["W_ad"] = r["m"] * r["w"]                     # potência de compressão adiabática
    pc = perda_compressor(t_sup(d), d["Tamb_C"], d["Patm_kPa"])
    r.update(h_conv=pc["h_conv"], h_rad=pc["h_rad"], h_comb=pc["h_conv"] + pc["h_rad"])
    r["Q_comp"] = pc["Q"] / 1e3                     # kW
    r["W_comp"] = r["W_ad"] + r["Q_comp"]           # potência de compressão
    r["COP"] = r["Q_evap"] / r["W_comp"]
    r["COP_ad"] = r["Q_evap"] / r["W_ad"]           # se o compressor fosse adiabático
    return r


def estados_ph(d):
    """(h, P) dos pontos medidos para o diagrama P-h."""
    P, r = pressoes(d), calcular(d)
    h = {2: r["h2"], 3: r["h3"], 6: r["h6"], 7: r["h7"], 8: h_TP(d["T8_C"], P[8]),
         9: h_Px(P[9], 0), 12: r["h12"]}                # ponto 9: título ~ 0
    return {k: (h[k], P[k]) for k in h}


# ---------------------------------------------------------------- incertezas
def _termos(funcao, d, u):
    """Contribuições (fator 1) de cada entrada, por diferenças centrais (GUM)."""
    termos = {}
    for nome, ui in u.items():
        dz = max(abs(d[nome]) * 1e-6, ui * 1e-3)
        dp, dm = dict(d), dict(d)
        dp[nome] += dz
        dm[nome] -= dz
        rp, rm = funcao(dp), funcao(dm)
        termos[nome] = {k: ((rp[k] - rm[k]) / (2 * dz) * ui) ** 2 for k in rp}
    return termos


def incerteza_numerica(funcao, d, u):
    """Incerteza expandida (fator de abrangência 2) de cada resultado."""
    termos = _termos(funcao, d, u)
    chaves = next(iter(termos.values())).keys()
    return {k: 2 * np.sqrt(sum(t[k] for t in termos.values())) for k in chaves}


def fracao_variancia(chave, entradas, d, u):
    """Parcela (%) da variância de `chave` que vem das `entradas`."""
    termos = _termos(calcular, d, u)
    total = sum(t[chave] for t in termos.values())
    return 100 * sum(termos[e][chave] for e in entradas) / total


# ---------------------------------------------------------------- conferência
def conferir_exemplo():
    """Exemplo do roteiro (R134a): o EES dá COP = 2,042 e W_comp = 0,3456 kW."""
    Patm, Pe, Pc = 100, 70 + 100, 1200 + 100
    hA = h_Px(Pc, 0, "R134a")
    hC, hD = h_TP(-5, Pe, "R134a"), h_TP(70, Pc, "R134a")
    m = 20 / 3600
    Qe, Qc = m * (hC - hA), 15 * 1200e-4 * (60 - 20) / 1000
    W = m * (hD - hC) + Qc
    print(f"Exemplo do roteiro (R134a): hA = {hA:.1f}, hC = {hC:.1f}, hD = {hD:.1f} kJ/kg; "
          f"Q_evap = {Qe:.4f} kW, W_comp = {W:.4f} kW, COP = {Qe / W:.3f}  "
          f"(EES: 122,6; 249,6; 298,9; 0,7058; 0,3456; 2,042)")


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


def _medida(v):
    """Valor medido com as casas da leitura (1 ou 2 decimais)."""
    return _br(v, 2 if round(v, 1) != v else 1)


def escrever_resultados_tex(d, r, U, arquivo="resultados.tex"):
    m = {"DataPratica": "03/11/2026 e 10/11/2026"}
    nomes = {"P2_man": "Pman", "P2_sat": "Psat", "x2": "Titulo", "SA3": "Superaq",
             "SR12": "Subresf", "h2": "Hdois", "h3": "Htres", "h6": "Hseis", "h7": "Hsete",
             "m": "Vazaokg", "q_evap": "EfeitoRef", "Q_evap": "CapRef", "w": "Wesp",
             "W_ad": "Wad", "h_conv": "Hconv", "h_rad": "Hrad", "h_comb": "Hcomb",
             "Q_comp": "Qcomp", "W_comp": "Wcomp", "COP": "COP", "COP_ad": "COPad"}
    escala = {"m": 1e3, "Q_comp": 1e3}             # g/s e W
    for k in CHAVES:
        a, b = _par(r[k] * escala.get(k, 1), U[k] * escala.get(k, 1))
        m[nomes[k]], m["U" + nomes[k]] = a, b

    # medições
    por_extenso = {1: "um", 2: "dois", 3: "tres", 4: "quatro", 5: "cinco", 6: "seis", 7: "sete",
                   8: "oito", 9: "nove", 10: "dez", 11: "onze", 12: "doze"}
    for k, s in por_extenso.items():
        m["T" + s] = _medida(d[f"T{k}_C"])
        if f"P{k}_kgf" in d:
            m["P" + s] = _medida(d[f"P{k}_kgf"])
    m.update({"Pdoze": _medida(d["P12_bar"]), "Vazao": _medida(d["m_lbm_min"]),
              "Tamb": _medida(d["Tamb_C"]), "Patm": _medida(d["Patm_kPa"]),
              "UPatm": _br(U_REL_PATM * d["Patm_kPa"], 2), "UPb": _br(U_PB, 2),
              "UPa": _br(U_PA, 2), "UPabar": _br(U_PA_BAR, 2), "UVazao": _br(U_M, 1)})

    # estados do ciclo (valores nominais)
    for k, (h, P) in estados_ph(d).items():
        m["Pabs" + por_extenso[k]] = _br(P, 0)
        m["Hp" + por_extenso[k]] = _br(h, 1)
    P = pressoes(d)
    pc = perda_compressor(t_sup(d), d["Tamb_C"], d["Patm_kPa"])
    ra = f"{pc['Ra']:.1e}".split("e")
    m.update({"Tsup": _br(pc["Tsup"], 1), "Sext": _br(pc["S"], 3),
              "Scm": _br(1e4 * pc["S"], 0), "Lcubo": _br(pc["L"], 3),
              "Rayleigh": f"{_br(float(ra[0]), 1)}\\times10^{{{int(ra[1])}}}",
              "TsatTres": _br(T_sat(P[3]), 1), "TsatDoze": _br(T_sat(P[12]), 1),
              "TsatNove": _br(T_sat(P[9]), 1), "TsatOito": _br(T_sat(P[8]), 1),
              "DifP": _br(abs(r["P2_man"] - r["P2_sat"]), 1),
              "FracQcomp": _br(100 * r["Q_comp"] / r["W_comp"], 0),
              "QcompW": _br(1e3 * r["Q_comp"], 0), "HconvNom": _br(r["h_conv"], 1),
              "HradNom": _br(r["h_rad"], 1), "HcombNom": _br(r["h_comb"], 1),
              "DPdT": _br((P_sat(d["T2_C"] + 0.5) - P_sat(d["T2_C"] - 0.5)) * U_T, 0),
              "UPbkPa": _br(U_PB * KGF, 0)})

    # faixa do COP conforme a temperatura da carcaça (de T_amb, sem perda, até T7)
    r_quente = calcular(dict(d, Tsup_C=d["T7_C"]))
    cop_c = (d["T2_C"] + 273.15) / (d["T9_C"] - d["T2_C"])
    u = incertezas(d)
    m.update({"COPquente": _br(r_quente["COP"], 1), "COPfrio": _br(r["COP_ad"], 1),
              "QcompQuente": _br(1e3 * r_quente["Q_comp"], 0),
              "COPcarnot": _br(cop_c, 1), "FracCarnot": _br(100 * r["COP"] / cop_c, 0),
              "FracVarW": _br(fracao_variancia(
                  "COP", ["T6_C", "T7_C", "P6_kgf", "P7_kgf"], d, u), 0),
              "FracVarM": _br(fracao_variancia("Q_evap", ["m_lbm_min"], d, u), 0)})

    with open(arquivo, "w", encoding="utf-8") as arq:
        arq.write("% Gerado automaticamente por codigo-07.py. Não editar à mão.\n")
        for nome, valor in m.items():
            arq.write(f"\\newcommand{{\\{nome}}}{{{valor}}}\n")
    print(f"\n{arquivo} gerado com {len(m)} macros.")


# ---------------------------------------------------------------- gráficos
def _txt(v, U):
    """'307 ± 47' com vírgula decimal, para os rótulos dos gráficos."""
    a, b = _par(v, U)
    return f"{a} ± {b}".replace("{,}", ",").replace("\\,", ".")


def grafico_ph(d, arquivo="ciclo_ph.png"):
    pontos = estados_ph(d)
    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    # curva de saturação
    Tc = CP.PropsSI("Tcrit", FLUIDO)
    T = np.linspace(233.15 - 10, Tc - 0.02, 400)
    hl = np.array([CP.PropsSI("H", "T", t, "Q", 0, FLUIDO) for t in T]) / 1e3
    hv = np.array([CP.PropsSI("H", "T", t, "Q", 1, FLUIDO) for t in T]) / 1e3
    ps = np.array([CP.PropsSI("P", "T", t, "Q", 0, FLUIDO) for t in T]) / 1e3
    ax.plot(np.r_[hl, hv[::-1]], np.r_[ps, ps[::-1]], color="0.15", lw=1.3,
            label="Curva de saturação do R12")

    # isotermas
    for i, Tiso in enumerate((0, 20, 40, 60, 80)):
        psat = P_sat(Tiso)
        pv = np.geomspace(100, psat * 0.999, 60)
        pl = np.geomspace(psat * 1.001, 4000, 30)
        h = np.r_[[h_TP(Tiso, p) for p in pv], h_Px(psat, 1), h_Px(psat, 0),
                  [h_TP(Tiso, p) for p in pl]]
        p = np.r_[pv, psat, psat, pl]
        ax.plot(h, p, color="0.72", lw=0.8, zorder=1,
                label="Isotermas (°C)" if i == 0 else None)
        ax.text(h_TP(Tiso, 160), 160, f"{Tiso}", color="0.45", fontsize=8,
                ha="center", va="center", bbox=dict(fc="white", ec="none", pad=0.6))

    # ciclo medido
    ordem = [12, 2, 3, 6, 7, 8, 9, 12]
    ax.plot([pontos[k][0] for k in ordem], [pontos[k][1] for k in ordem], "-o",
            color="C0", lw=2, ms=6, mec="white", mew=1.2, zorder=3,
            label="Ciclo da bancada (pontos medidos)")
    desloc = {12: (-14, 6), 2: (-12, -14), 3: (-6, -16), 6: (8, -6), 7: (8, 4),
              8: (-12, 8), 9: (4, 8)}
    for k, (h, P) in pontos.items():
        ax.annotate(str(k), (h, P), textcoords="offset points", xytext=desloc[k],
                    fontsize=10, fontweight="bold", color="0.1")
    rotulo = dict(fontsize=8, color="0.35", style="italic",
                  bbox=dict(fc="white", ec="none", pad=0.8, alpha=0.9))
    ax.text(pontos[12][0] + 3, 430, "válvula\ntermostática", ha="left", va="center", **rotulo)
    ax.text((pontos[2][0] + pontos[3][0]) / 2, pontos[2][1] * 0.9, "evaporador",
            ha="center", va="top", **rotulo)
    ax.text((pontos[8][0] + pontos[9][0]) / 2, pontos[9][1] * 1.08, "condensador",
            ha="center", va="bottom", **rotulo)
    p_meio = np.sqrt(pontos[6][1] * pontos[7][1])
    h_meio = np.interp(np.log(p_meio), np.log([pontos[6][1], pontos[7][1]]),
                       [pontos[6][0], pontos[7][0]])
    ax.text(h_meio + 4, p_meio, "compressor", ha="left", va="center", **rotulo)

    ax.set_yscale("log")
    ax.set_xlim(0, 260)
    ax.set_ylim(150, 3000)
    ax.yaxis.set_major_locator(FixedLocator([200, 300, 500, 1000, 2000, 3000]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
    ax.set_xlabel("Entalpia específica, $h$ (kJ/kg)")
    ax.set_ylabel("Pressão absoluta, $P$ (kPa)")
    ax.grid(alpha=0.25, which="major")
    ax.legend(fontsize=8, loc="upper left", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(arquivo, dpi=300)
    return fig


def grafico_pressao(r, U, arquivo="pressao_entrada.png"):
    linhas = [("Saturação a $T_2$", r["P2_sat"], U["P2_sat"], "C1"),
              ("Manômetro + $P_{atm}$", r["P2_man"], U["P2_man"], "C0")]
    fig, ax = plt.subplots(figsize=(7.6, 2.3))
    lo = max(v - u for _, v, u, _ in linhas)
    hi = min(v + u for _, v, u, _ in linhas)
    ax.axvspan(lo, hi, color="0.9", zorder=0)
    ax.text((lo + hi) / 2, -0.62, "interseção das faixas", ha="center", va="center",
            fontsize=8, color="0.35")
    for y, (nome, v, u, cor) in enumerate(linhas):
        ax.barh(y, 2 * u, left=v - u, height=0.36, color=cor, alpha=0.28, edgecolor=cor, lw=1)
        ax.errorbar(v, y, xerr=u, fmt="o", color=cor, ms=7, capsize=5, lw=1.5)
        ax.text(v, y + 0.26, f"{_txt(v, u)} kPa", ha="center", va="bottom", fontsize=9)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([n for n, *_ in linhas])
    ax.set_ylim(-0.85, 1.75)
    ax.set_xlim(min(v - u for _, v, u, _ in linhas) - 8, max(v + u for _, v, u, _ in linhas) + 8)
    ax.set_xlabel("Pressão absoluta na entrada do evaporador (kPa)")
    ax.grid(alpha=0.25, axis="x")
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    fig.tight_layout()
    fig.savefig(arquivo, dpi=300)
    return fig


# ---------------------------------------------------------------- principal
def main():
    conferir_exemplo()
    d = DADOS
    r = calcular(d)
    U = incerteza_numerica(calcular, d, incertezas(d))

    with open("resultados_pratica7.csv", "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(["grandeza", "valor", "U (fator 2)"])
        w.writerows((k, round(r[k], 6), round(U[k], 6)) for k in CHAVES)
    print()
    for k in CHAVES:
        print(f"{k:8s} = {r[k]:12.5g} ± {U[k]:.2g}")
    for k, (h, P) in estados_ph(d).items():
        print(f"ponto {k:2d}: P = {P:7.1f} kPa, h = {h:6.1f} kJ/kg")

    escrever_resultados_tex(d, r, U)
    grafico_pressao(r, U)
    grafico_ph(d)
    plt.show()


if __name__ == "__main__":
    main()
