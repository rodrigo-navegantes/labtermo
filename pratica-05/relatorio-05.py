\documentclass[11pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[brazil]{babel}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{geometry}
\geometry{top=2.0cm,bottom=2.0cm,left=2.2cm,right=2.2cm,headheight=14pt}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs,tabularx,array,multirow}
\usepackage{xcolor}
\usepackage{fancyhdr}
\usepackage{lastpage}
\usepackage{caption}
\usepackage{float}
\usepackage{hyperref}

\definecolor{petrol}{HTML}{1F4E5F}

\hypersetup{
  colorlinks=true,
  linkcolor=petrol,
  citecolor=petrol,
  urlcolor=petrol,
  pdftitle={Transferência de calor por ebulição em vaso},
  pdfauthor={Rodrigo Malato Navegantes}
}

\pagestyle{fancy}
\fancyhf{}
\lhead{EMA-103 -- Laboratório de Térmica}
\rhead{Prática 5}
\cfoot{\thepage\ de \pageref{LastPage}}
\renewcommand{\headrulewidth}{0.4pt}

\setlength{\parindent}{1.2em}
\setlength{\parskip}{0.15em}
\renewcommand{\arraystretch}{1.12}
\captionsetup{font=small,labelfont=bf}
\renewcommand{\refname}{Referências}

\makeatletter
\renewcommand\section{\@startsection{section}{1}{\z@}%
  {-2.4ex \@plus -1ex \@minus -.2ex}{1.0ex \@plus .2ex}%
  {\normalfont\large\bfseries\color{petrol}}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
  {-1.8ex \@plus -1ex \@minus -.2ex}{0.6ex \@plus .2ex}%
  {\normalfont\normalsize\bfseries\color{petrol}}}
\makeatother

\input{resultados.tex}

\begin{document}

\begin{center}
  {\small UNIVERSIDADE FEDERAL DE MINAS GERAIS\\
  Escola de Engenharia -- Engenharia Aeroespacial\\
  EMA-103 -- Laboratório de Térmica}

  \vspace{1.25cm}
  {\LARGE\bfseries\color{petrol} Transferência de calor por ebulição em vaso}

  \vspace{0.95cm}
  \begin{tabular}{>{\bfseries}r l}
    Aluno(a): & Rodrigo Malato Navegantes\\
    Matrícula: & 2023038175\\
    Turma: & EA1\\
    Professor(a): & Luiz Machado\\
    Data da prática: & \DataPratica\\
  \end{tabular}
\end{center}

\vspace{0.55cm}\hrule\vspace{0.65cm}

\section{Objetivo}

Determinar o coeficiente de transferência de calor por ebulição em vaso,
$h_{eb}$, a partir de medições em um calefator de cobre imerso em água fervente,
e comparar o valor obtido com as correlações de Rohsenow e de Foster/Zuber.

\section{Materiais e montagem}

Foram utilizados: calefator de cobre aquecido por resistores elétricos e
instrumentado com três termopares; vaso de pressão com água fervente,
manômetro e termopar; condensador de vidro resfriado com água da rede;
termopares tipo K e medidor de temperatura; multímetro (tensão $V_{AB}$);
amperímetro alicate (corrente $I$); escala de aço; paquímetro; bureta;
cronômetro. O esquema e a foto da montagem estão nas
Figuras~\ref{fig:esquema} e~\ref{fig:foto}.

\begin{figure}[H]
  \centering
  \includegraphics[width=0.86\textwidth]{esquema5.png}
  \caption{Esquema da montagem: o calor gerado pelos resistores atravessa o
  calefator de cobre (termopares $T_1$, $T_2$ e $T_3$) e ferve a água do vaso; o
  vapor é condensado no condensador de vidro. Fonte: material didático da
  disciplina~\cite{roteiro}.}
  \label{fig:esquema}
\end{figure}

\begin{figure}[H]
  \centering
  \includegraphics[width=0.74\textwidth]{foto5.jpg}
  \caption{Foto da montagem. Fonte: material didático da
  disciplina~\cite{roteiro}.}
  \label{fig:foto}
\end{figure}

\section{Metodologia}

\subsection{Coeficiente experimental}

Em regime permanente, a potência elétrica dos resistores atravessa o topo do
calefator, de área $S=\pi d^2/4$, e é transferida à água por ebulição:
\begin{equation}
  q=V_{AB}\,I=h_{eb}\,S\,(T_s-T_{eb})
  \quad\Rightarrow\quad
  h_{eb}=\frac{V_{AB}\,I}{S\,\Delta T_e},\qquad \Delta T_e=T_s-T_{eb}.
  \label{eq:hexp}
\end{equation}
A temperatura da superfície $T_s$ não é medida diretamente. Como o calefator é
isolado lateralmente, a condução no cobre é unidimensional e a temperatura varia
linearmente com a profundidade $x$; por isso, ajustou-se uma reta às leituras
$T_1$, $T_2$ e $T_3$, nas profundidades $x_1$, $x_2$ e $x_3$, e extrapolou-se
até o topo ($x=0$).

\subsection{Correlações}

Pela correlação de Rohsenow~\cite{roteiro,incropera},
\begin{equation}
  h_{eb}=\mu_l\,h_{lv}\left[\frac{g(\rho_l-\rho_v)}{\sigma}\right]^{1/2}
  \left(\frac{c_{p,l}\,\Delta T_e^{2/3}}{C_{s,f}\,h_{lv}\,\mathrm{Pr}_l^{\,n}}\right)^3,
  \label{eq:rohsenow}
\end{equation}
com $n=1$ e, para água e cobre, $C_{s,f}=0{,}0128$ (cobre polido) ou
$C_{s,f}=0{,}0068$ (cobre riscado). Como o estado da superfície do calefator não
é conhecido, os dois valores foram avaliados. A correlação de Foster/Zuber foi
aplicada na forma do roteiro,
\begin{equation}
  h_{eb}=1{,}22\left(\frac{k_l^{0,79}\,c_{p,l}^{0,45}\,\rho_l^{0,49}}
  {\sigma^{0,5}\,\mu_l^{0,29}\,h_{lv}^{0,24}\,\rho_v^{0,24}}\right)^3
  \Delta T_e^{0,24}\,\Delta P_{sat}^{0,75},
  \label{eq:fz}
\end{equation}
com $k_l$ em kW/(m\,K), $c_{p,l}$ em kJ/(kg\,K), $h_{lv}$ em kJ/kg e
$\Delta P_{sat}=P_{sat}(T_s)-P_{sat}(T_{eb})$ em Pa. As propriedades da água
saturada foram avaliadas em $T_{eb}$ (Tabela~\ref{tab:props}) com a biblioteca
CoolProp, que usa a formulação IAPWS.

As incertezas foram propagadas numericamente pelo método do GUM, a partir das
incertezas da tabela de medições do roteiro, e são apresentadas como incertezas
expandidas (fator de abrangência 2).

\section{Resultados e discussão}

\begin{table}[H]
\centering
\caption{Medições. A pressão atmosférica é a da atmosfera padrão na altitude
de Belo Horizonte (cerca de 860 m).}
\label{tab:medicoes}
\small
\begin{tabular}{lcc@{\hspace{2.2em}}lcc}
\toprule
\textbf{Grandeza} & \textbf{Valor} & \textbf{Incerteza} &
\textbf{Grandeza} & \textbf{Valor} & \textbf{Incerteza}\\
\midrule
$d$ (mm)      & 30   & $\pm0{,}1$ mm & $T_{eb}$ ($^{\circ}$C) & \Teb & $\pm1\,^{\circ}$C\\
$x_1$ (mm)    & 3    & $\pm0{,}5$ mm & $V_{cond}$ (mL) & \Vcond & $\pm1\%$\\
$x_2$ (mm)    & 22   & $\pm0{,}1$ mm & $\Delta t_{cond}$ (s) & \tcond & $\pm5$ s\\
$x_3$ (mm)    & 31   & $\pm0{,}1$ mm & $T_{a1}$ ($^{\circ}$C) & \Taum & $\pm1\,^{\circ}$C\\
$T_1$ ($^{\circ}$C) & \Tum & $\pm1\,^{\circ}$C & $T_{a2}$ ($^{\circ}$C) & \Tadois & $\pm1\,^{\circ}$C\\
$T_2$ ($^{\circ}$C) & \Tdois & $\pm1\,^{\circ}$C & $V_{resf}$ (mL) & \Vresf & $\pm1\%$\\
$T_3$ ($^{\circ}$C) & \Ttres & $\pm1\,^{\circ}$C & $\Delta t_{resf}$ (s) & \tresf & $\pm5$ s\\
$V_{AB}$ (V)  & \VAB & $\pm5\%$ & $T_{\mathrm{amb}}$ ($^{\circ}$C) & \Tamb & $\pm1\,^{\circ}$C\\
$I$ (A)       & \Corrente & $\pm5\%$ & $P_{atm}$ (kPa) & \Patm & $\pm1\%$\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[H]
\centering
\caption{Propriedades da água saturada em $T_{eb}=\Teb\,^{\circ}$C.}
\label{tab:props}
\small
\begin{tabular}{cccccccc}
\toprule
$\rho_l$ & $\rho_v$ & $\mu_l$ & $k_l$ & $c_{p,l}$ & $\mathrm{Pr}_l$ & $\sigma$ & $h_{lv}$\\
(kg/m$^3$) & (kg/m$^3$) & ($10^{-4}$ Pa\,s) & (W/(m\,K)) & (J/(kg\,K)) & -- & (mN/m) & (kJ/kg)\\
\midrule
\Rhol & \Rhov & \Mul & \Kl & \Cpl & \Prl & \Tensao & \Hlv\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[H]
\centering
\caption{Resultados (incertezas expandidas, fator de abrangência 2).}
\label{tab:resultados}
\small
\begin{tabular}{lcc}
\toprule
\textbf{Grandeza} & \textbf{Valor} & \textbf{Razão ao experimental}\\
\midrule
Temperatura da superfície, $T_s$ ($^{\circ}$C) & $\Ts\pm\UTs$ & \\
Excesso de temperatura, $\Delta T_e$ (K) & $\DTe\pm\UDTe$ & \\
Potência elétrica, $q=V_{AB}I$ (W) & $\Qel\pm\UQel$ & \\
Fluxo de calor, $q''=q/S$ (kW/m$^2$) & $\Fluxo\pm\UFluxo$ & \\
\midrule
$h_{eb}$ experimental (W/(m$^2$\,K)) & $\Hexp\pm\UHexp$ & 1\\
$h_{eb}$ Rohsenow, cobre polido (W/(m$^2$\,K)) & $\Hpolido\pm\UHpolido$ & \RazaoPolido\\
$h_{eb}$ Rohsenow, cobre riscado (W/(m$^2$\,K)) & $\Hriscado\pm\UHriscado$ & \RazaoRiscado\\
$h_{eb}$ Foster/Zuber (W/(m$^2$\,K)) & $\Hfz\pm\UHfz$ & \RazaoFz\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}[H]
  \centering
  \IfFileExists{ebulicao.png}{\includegraphics[width=\textwidth]{ebulicao.png}}%
    {\fbox{\parbox[c][4cm][c]{0.9\textwidth}{\centering
     Rode o codigo-05.py com as medições para gerar ebulicao.png.}}}
  \caption{(a) Temperaturas medidas no calefator e extrapolação linear até o
  topo ($x=0$), que fornece $T_s$. (b) Fluxo de calor em função do excesso de
  temperatura: correlações e ponto experimental, com barras de incerteza
  expandida.}
  \label{fig:ebulicao}
\end{figure}

A temperatura de ebulição medida, $\Teb\,^{\circ}$C, é compatível com a de
saturação na pressão atmosférica local, $\Tsatatm\,^{\circ}$C. Com
$\Delta T_e=\DTe\pm\UDTe$ K, a ebulição ocorreu no regime nucleado, na região de
jatos e colunas da curva de ebulição ($10<\Delta T_e<30$ K)~\cite{incropera}.

O coeficiente experimental, $h_{eb}=\Hexp\pm\UHexp$ W/(m$^2$\,K), ficou abaixo
de todas as correlações (Tabela~\ref{tab:resultados} e
Figura~\ref{fig:ebulicao}b): Foster/Zuber o superestima por um fator de
\RazaoFz, e Rohsenow por \RazaoPolido{} (cobre polido) e \RazaoRiscado{} (cobre
riscado). A diferença é coerente com a sensibilidade dessas correlações: em
Rohsenow, $h_{eb}$ varia com $C_{s,f}^{-3}$ e com $\Delta T_e^{2}$, de modo que
pequenas diferenças no estado da superfície ou em $T_s$ mudam muito a previsão,
e erros de até 100\% no fluxo de calor são esperados~\cite{incropera}. O valor
de $C_{s,f}$ que reproduziria a medição é \CsfAjuste, maior que o do cobre
polido, o que indica uma superfície com menos pontos de nucleação que a de
referência.

Nem toda a potência elétrica chega à água em ebulição. O condensado coletado
corresponde a $\dot m_{cond}\,h_{lv}=\Qevap\pm\UQevap$ W, cerca de
\Fracevap\% de $V_{AB}I$, ordem de grandeza confirmada pelo calor retirado pela
água de resfriamento ($\Qresf\pm\UQresf$ W) e pela condução no próprio cobre,
$k\,(dT/dx)\,S=\Qcond\pm\UQcond$ W (com $k=395$ W/(m\,K)). O restante se perde
pelo isolamento e pelas paredes do vaso. Assim, o $h_{eb}$ calculado com
$V_{AB}I$ é um limite superior: o valor real seria ainda menor, o que reforça
que, nesta montagem, as correlações superestimam o coeficiente de ebulição.

Todo o código empregado nesta prática esta disposto no \href{Repositório GIT da disciplina}{}

\begin{thebibliography}{9}
\bibitem{roteiro}
UFMG. \textit{EMA-103: Laboratório de Térmica -- Aulas 10 e 11: transferência
de calor por ebulição}. Material didático da disciplina, 2026.

\bibitem{incropera}
BERGMAN, T. L.; LAVINE, A. S.; INCROPERA, F. P.; DEWITT, D. P.
\textit{Fundamentals of Heat and Mass Transfer}. 7. ed. Hoboken: Wiley, 2011.
\end{thebibliography}

\end{document}
