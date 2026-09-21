"""Prática 3 — aleta cilíndrica. Preencha somente o bloco DADOS."""

import csv

import matplotlib.pyplot as plt
import numpy as np


# ================================ DADOS =====================================
DADOS = {
    "L_mm": 297.0,
    "d_mm": 25.4,
    "Tamb_C": 23.3,
    # Temperaturas medidas em x = 0, 45, 95, 145, 195, 245 e 295 mm.
    "T_C": [94.1, 79.0, 66.8, 55.9, 53.1, 49.8, 48.4],
}
# ============================================================================

X_MM = np.array([0, 45, 95, 145, 195, 245, 295], dtype=float)
K, U_K = 60.0, 5.0                  # W/(m K)
H, U_H = 12.0, 2.0                  # W/(m² K)
U_L_MM, U_D_MM, U_T_C = 0.5, 0.1, 1.0

SUFIXOS = [
    "zero", "quarentacinco", "noventacinco", "centoquarentacinco",
    "centonoventacinco", "duzentosquarentacinco", "duzentosnoventacinco",
]


def calcular(z, x=None):
    """Retorna geometria, perfil teórico, taxas, efetividade e eficiência."""
    L, d = z[:2]
    T = z[2:9]
    Tamb, k, h = z[9:]
    x = X_MM / 1000 if x is None else x

    P, A = np.pi * d, np.pi * d**2 / 4
    m = np.sqrt(h * P / (k * A))
    b = h / (m * k)
    denominador = np.cosh(m * L) + b * np.sinh(m * L)

    Tteo = Tamb + (T[0] - Tamb) * (
        np.cosh(m * (L - x)) + b * np.sinh(m * (L - x))
    ) / denominador
    qteo = np.sqrt(h * P * k * A) * (T[0] - Tamb) * (
        np.sinh(m * L) + b * np.cosh(m * L)
    ) / denominador

    # O último ponto está 2 mm antes da ponta: adota-se T(295 mm) até L.
    x_int = np.r_[X_MM / 1000, L]
    theta = np.r_[T, T[-1]] - Tamb
    qexp = h * P * np.trapz(theta, x_int) + h * A * theta[-1]

    Af = P * L + A
    theta_b = T[0] - Tamb
    efet_teo, efet_exp = qteo / (h * A * theta_b), qexp / (h * A * theta_b)
    efic_teo, efic_exp = qteo / (h * Af * theta_b), qexp / (h * Af * theta_b)
    return np.r_[P, A, m, Af, Tteo, qteo, qexp, efet_teo, efet_exp, efic_teo, efic_exp]


def incerteza(funcao, z, uz):
    """Propagação numérica por diferenças centrais; retorna U com k=2."""
    variancia = 0.0
    for i, ui in enumerate(uz):
        passo = max(abs(z[i]) * 1e-6, ui * 1e-3, 1e-8)
        maior, menor = z.copy(), z.copy()
        maior[i], menor[i] = maior[i] + passo, menor[i] - passo
        derivada = (funcao(maior) - funcao(menor)) / (2 * passo)
        variancia += (derivada * ui) ** 2
    return 2 * np.sqrt(variancia)


def br(x, casas=1):
    return f"{x:.{casas}f}".replace(".", "{,}")


def escrever_tex(z, y, U, arquivo="resultados.tex"):
    """Gera as macros usadas por relatorio_simplificado.tex."""
    n = len(X_MM)
    Tteo, UTteo = y[4:4 + n], U[4:4 + n]
    qteo, qexp, efet_teo, _, efic_teo, _ = y[4 + n:]
    Uqteo, Uqexp, Uefet_teo, _, Uefic_teo, _ = U[4 + n:]
    Texp = np.asarray(DADOS["T_C"], dtype=float)

    delta = Tteo - Texp
    compativel = np.abs(delta) <= UTteo + 2 * U_T_C
    criterio = np.abs(delta) <= 1.0
    maior = 1 + np.argmax(np.abs(delta[1:]))
    por_extenso = ["nenhum", "um", "dois", "três", "quatro", "cinco", "seis"]

    macros = {
        "Laleta": br(DADOS["L_mm"], 0),
        "Daleta": br(DADOS["d_mm"], 1),
        "Tamb": br(DADOS["Tamb_C"], 1),
        "mL": br(y[2] * z[0], 3),
        "Qteorico": br(qteo), "UQteorico": br(Uqteo),
        "Qexperimental": br(qexp), "UQexperimental": br(Uqexp),
        "DesvioQ": br(100 * (qexp - qteo) / qteo),
        "EfetividadeTeorica": br(efet_teo),
        "UEfetividadeTeorica": br(Uefet_teo),
        "EficienciaTeorica": br(efic_teo, 3),
        "UEficienciaTeorica": br(Uefic_teo, 3),
        "NCompativeis": por_extenso[int(compativel[1:].sum())],
        "NRoteiro": por_extenso[int(criterio[1:].sum())],
        "XMaiorDesvio": f"{X_MM[maior]:g}",
        "MaiorDesvio": br(abs(delta[maior])),
        "QualificacaoPerfil": "boa" if compativel[1:].sum() >= 5 else "apenas parcial",
        "QualificacaoTaxa": (
            "compatíveis" if abs(qexp - qteo) <= Uqteo + Uqexp else "incompatíveis"
        ),
    }

    for i, sufixo in enumerate(SUFIXOS):
        macros[f"T{sufixo}"] = br(Texp[i])
        macros[f"Tteo{sufixo}"] = br(Tteo[i])
        macros[f"UTteo{sufixo}"] = br(UTteo[i])
        macros[f"Delta{sufixo}"] = "--" if i == 0 else (
            "$" + ("+" if delta[i] > 0 else "-") + br(abs(delta[i])) + "$"
        )
        macros[f"Comp{sufixo}"] = "--" if i == 0 else (
            "Sim" if compativel[i] else "Não"
        )

    with open(arquivo, "w", encoding="utf-8") as arq:
        arq.write("% Gerado automaticamente por codigo-03.py. Não editar à mão.\n")
        for nome, valor in macros.items():
            arq.write(f"\\newcommand{{\\{nome}}}{{{valor}}}\n")


def main():
    T = np.asarray(DADOS["T_C"], dtype=float)
    if len(T) != len(X_MM) or not np.isfinite(T).all():
        raise ValueError("Informe as sete temperaturas medidas em T_C.")
    if DADOS["L_mm"] <= X_MM[-1] or DADOS["d_mm"] <= 0:
        raise ValueError("Use L_mm > 295 mm e d_mm positivo.")

    z = np.r_[DADOS["L_mm"] / 1000, DADOS["d_mm"] / 1000, T, DADOS["Tamb_C"], K, H]
    uz = np.r_[U_L_MM / 1000, U_D_MM / 1000, [U_T_C] * 8, U_K, U_H]
    y, U = calcular(z), incerteza(calcular, z, uz)

    nomes = (
        ["P (m)", "A (m²)", "m (1/m)", "Af (m²)"]
        + [f"Tteo({x:g} mm) (°C)" for x in X_MM]
        + ["qteo (W)", "qexp (W)", "efet_teo", "efet_exp", "efic_teo", "efic_exp"]
    )
    with open("resultados_pratica3.csv", "w", newline="", encoding="utf-8") as arq:
        escritor = csv.writer(arq)
        escritor.writerow(["grandeza", "valor", "U (k=2)"])
        escritor.writerows(zip(nomes, y, U))
    escrever_tex(z, y, U)

    x = np.linspace(0, z[0], 300)
    perfil = lambda zz: calcular(zz, x)[4:4 + len(x)]
    Tteo, UTteo = perfil(z), incerteza(perfil, z, uz)
    plt.figure(figsize=(6.4, 4.2))
    plt.plot(1000 * x, Tteo, color="C1", label="Modelo analítico")
    plt.fill_between(1000 * x, Tteo - UTteo, Tteo + UTteo, color="C1", alpha=.18,
                     lw=0, label="Modelo, $U$ com $k=2$")
    plt.errorbar(X_MM, T, yerr=2 * U_T_C, fmt="o", color="C0", capsize=3,
                 zorder=3, label="Experimental ($k=2$)")
    plt.xlabel("Posição, $x$ (mm)")
    plt.ylabel("Temperatura (°C)")
    plt.grid(alpha=.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig("perfil_temperatura.png", dpi=300)
    plt.close()

    print("Gerados: resultados_pratica3.csv, resultados.tex e perfil_temperatura.png.")


if __name__ == "__main__":
    main()

