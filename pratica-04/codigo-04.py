"""Prática 4 - aleta cônica"""

import csv
import numpy as np
import matplotlib.pyplot as plt

DADOS = {
    "Tb_C": 95.3,       # T(x = 0) 
    "T60_C": 80.4,      # T(x = 60 mm)
    "T120_C": 69.0,     # T(x = 120 mm)
    "T180_C": 58.4,     # T(x = 180 mm)
    "T240_C": 49.9,     # T(x = 240 mm)
    "Tamb_C": 27.5,
}
X_MM = np.array([0, 60, 120, 180, 240], dtype=float)


L_MM, D0_MM, DL_MM = 297.0, 25.4, 1.5
U_L_MM, U_D0_MM, U_DL_MM, U_T_C = 0.5, 0.1, 0.1, 1.0


K, H = 60.0, 12.0          
U_K, U_H = 5.0, 2.0
N = 20                   


def perfil_mdf(L, d0, dL, Tb, Tamb, k=K, h=H, n=N):

    dx = L / n
    x = np.arange(n + 1) * dx
    raio = lambda xx: d0 / 2 + (dL - d0) / 2 * xx / L
    r1 = raio(x[1:] - dx / 2)                  # face voltada para a base
    r2 = raio(np.minimum(x[1:] + dx / 2, L))   # face voltada para a ponta
    A1, A2 = np.pi * r1 ** 2, np.pi * r2 ** 2
    Sx = np.pi * (r1 + r2) * dx
    c1, c2 = k * A1 / dx, k * A2 / dx
    A_ponta = np.pi * (dL / 2) ** 2

    M, b = np.zeros((n, n)), np.zeros(n)
    for j in range(n):                         # linha j <-> nó i = j + 1
        ponta = j == n - 1
        perda = h * (Sx[j] / 2 + A_ponta) if ponta else h * Sx[j]
        M[j, j] = c1[j] + perda + (0.0 if ponta else c2[j])
        if j > 0:
            M[j, j - 1] = -c1[j]
        if not ponta:
            M[j, j + 1] = -c2[j]
        b[j] = perda * Tamb + (c1[j] * Tb if j == 0 else 0.0)
    return 1000 * x, np.r_[Tb, np.linalg.solve(M, b)]


def calcular(z):
    x_mm, T = perfil_mdf(*z)
    return np.interp(X_MM, x_mm, T)


def incerteza_numerica(funcao, valores, incertezas, separar=False):
    termos = []
    for i, ui in enumerate(incertezas):
        dz = max(abs(valores[i]) * 1e-6, ui * 1e-3, 1e-8)
        zp, zm = valores.copy(), valores.copy()
        zp[i], zm[i] = zp[i] + dz, zm[i] - dz
        termos.append(((funcao(zp) - funcao(zm)) * ui / (2 * dz)) ** 2)
    termos = np.array(termos)
    U = 2 * np.sqrt(termos.sum(axis=0))
    return (U, termos) if separar else U


def h_ajuste(z, Texp):
    hs = np.linspace(5, 30, 2501)
    sse = [np.sum((calcular(np.r_[z[:6], h]) - Texp)[1:] ** 2) for h in hs]
    return hs[int(np.argmin(sse))]


SUFIXOS = ["zero", "sessenta", "centovinte", "centoitenta", "duzentosquarenta"]


def _br(valor, casas=1):
    """Formata um número no padrão brasileiro, pronto para o LaTeX."""
    return f"{valor:.{casas}f}".replace(".", "{,}")


def escrever_resultados_tex(z, Texp, Tnum, Unum, termos, arquivo="resultados.tex"):
    """Gera o arquivo de macros consumido por Relatorio_Pratica_4.tex."""
    delta = Tnum - Texp
    compat = np.abs(delta) <= Unum + 2 * U_T_C
    ncomp = int(compat[1:].sum())          # x = 0 é condição de contorno
    imaior = 1 + int(np.argmax(np.abs(delta[1:])))
    hb = h_ajuste(z, Texp)
    desvio_aj = np.abs(calcular(np.r_[z[:6], hb]) - Texp)[1:].max()
    fracao_h = 100 * termos[6, -1] / termos[:, -1].sum()   # em x = 240 mm

    macros = {
        "DataPratica": "22/09/2026 e 29/09/2026",
        "Laleta": _br(L_MM, 0), "DBasealeta": _br(D0_MM, 1), "DLaleta": _br(DL_MM, 1),
        "Tamb": _br(DADOS["Tamb_C"], 1), "N": f"{N}",
        "Biot": _br(1e3 * H * (D0_MM / 2000) / K, 1),
        "FraseCompat": {0: "nenhum dos quatro pontos independentes tem",
                        1: "um dos quatro pontos independentes tem",
                        2: "dois dos quatro pontos independentes têm",
                        3: "três dos quatro pontos independentes têm",
                        4: "os quatro pontos independentes têm"}[ncomp],
        "XMaiorDesvio": f"{X_MM[imaior]:g}",
        "MaiorDesvio": _br(abs(delta[imaior]), 1),
        "FracaoH": _br(fracao_h, 0),
        "HAjuste": _br(hb, 1),
        "DesvioAjuste": _br(desvio_aj, 1),
    }
    for i, suf in enumerate(SUFIXOS):
        macros[f"Texp{suf}"] = _br(Texp[i], 1)
        macros[f"Tnum{suf}"] = _br(Tnum[i], 1)
        macros[f"UTnum{suf}"] = _br(Unum[i], 1)
        macros[f"Delta{suf}"] = ("--" if i == 0 else
                                  f"${'+' if delta[i] > 0 else '-'}{_br(abs(delta[i]), 1)}$")
        macros[f"Comp{suf}"] = "--" if i == 0 else ("Sim" if compat[i] else "Não")

    with open(arquivo, "w", encoding="utf-8") as arq:
        arq.write("% Gerado automaticamente por codigo-04.py. Não editar à mão.\n")
        for nome, valor in macros.items():
            arq.write(f"\\newcommand{{\\{nome}}}{{{valor}}}\n")
    print(f"\n{arquivo} gerado com {len(macros)} macros "
          f"(h de melhor ajuste = {hb:.1f} W/(m² K)).")


def main():
    Texp = np.array([DADOS["Tb_C"], DADOS["T60_C"], DADOS["T120_C"],
                     DADOS["T180_C"], DADOS["T240_C"]])
    if np.isnan(np.r_[Texp, DADOS["Tamb_C"]]).any():
        raise SystemExit("Preencha Tb_C, T60_C, T120_C, T180_C, T240_C e Tamb_C em DADOS.")

    z = np.r_[L_MM / 1000, D0_MM / 1000, DL_MM / 1000, DADOS["Tb_C"], DADOS["Tamb_C"], K, H]
    uz = np.r_[U_L_MM / 1000, U_D0_MM / 1000, U_DL_MM / 1000, U_T_C, U_T_C, U_K, U_H]
    Tnum = calcular(z)
    Unum, termos = incerteza_numerica(calcular, z, uz, separar=True)

    with open("resultados_pratica4.csv", "w", newline="", encoding="utf-8") as arq:
        w = csv.writer(arq)
        w.writerow(["x (mm)", "T_experimental (°C)", "T_MDF (°C)", "U_MDF (k=2)"])
        w.writerows(zip(X_MM, Texp, np.round(Tnum, 2), np.round(Unum, 2)))
    for x, te, tn, u in zip(X_MM, Texp, Tnum, Unum):
        print(f"x = {x:5.0f} mm   Texp = {te:5.1f} °C   T_MDF = {tn:5.1f} ± {u:3.1f} °C")

    escrever_resultados_tex(z, Texp, Tnum, Unum, termos)

    # Gráfico: perfil do MDF (nós), faixa de incerteza e pontos medidos.
    x_mm, Tcurva = perfil_mdf(*z)
    Ucurva = incerteza_numerica(lambda zz: perfil_mdf(*zz)[1], z, uz)
    plt.figure(figsize=(6.4, 4.2))
    plt.plot(x_mm, Tcurva, "-o", color="C1", ms=3, lw=1.6, zorder=2, label=f"MDF ($N={N}$)")
    plt.fill_between(x_mm, Tcurva - Ucurva, Tcurva + Ucurva, color="C1",
                     alpha=.18, lw=0, zorder=1, label="MDF, incerteza expandida")
    plt.errorbar(X_MM, Texp, yerr=2 * U_T_C, fmt="o", color="C0", ms=6,
                 capsize=3, zorder=3, label="Experimental")
    plt.xlabel("Posição, $x$ (mm)")
    plt.ylabel("Temperatura (°C)")
    plt.grid(alpha=.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("perfil_temperatura_conica.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
