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
  pdftitle={Transferência de calor por condensação},
  pdfauthor={Rodrigo Malato Navegantes}
}

\pagestyle{fancy}
\fancyhf{}
\lhead{EMA-103 -- Laboratório de Térmica}
\rhead{Prática 6}
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
  {\LARGE\bfseries\color{petrol} Transferência de calor por condensação}

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

Determinar o coeficiente de transferência de calor por condensação interna,
$h_{cond}$, na serpentina de um condensador de vidro, a partir de medições, e
comparar o resultado com a correlação de Dobson e Chato.

\section{Materiais e montagem}

A montagem é a mesma da prática de ebulição (Figura~\ref{fig:esquema}): a água
ferve no vaso de pressão, aquecida pelo calefator de cobre, e o vapor segue para
o condensador de vidro. Nele, o vapor condensa dentro de uma serpentina
helicoidal, resfriada pela água da rede que escoa no tubo envelope
(Figura~\ref{fig:condensador}). Foram utilizados: calefator de cobre com
resistores elétricos; vaso de pressão com água fervente; condensador de vidro
resfriado com água da rede; termopares tipo K e medidor de temperatura;
multímetro; amperímetro alicate; escala de aço; paquímetro; bureta; cronômetro.

\begin{figure}[H]
  \centering
  \includegraphics[width=0.84\textwidth]{esquema6.png}
  \caption{Esquema da montagem. Fonte: material didático da
  disciplina~\cite{roteiro}.}
  \label{fig:esquema}
\end{figure}

\begin{figure}[H]
  \centering
  \includegraphics[width=0.92\textwidth]{condensador6.jpg}
  \caption{Condensador: tubo envelope ($t$), hélice ($h$) e tubo da serpentina
  ($s$), com as dimensões medidas. O vapor condensa dentro da serpentina, e a
  água de resfriamento escoa no tubo envelope. Fonte: material didático da
  disciplina~\cite{roteiro}.}
  \label{fig:condensador}
\end{figure}

\section{Metodologia}

O cálculo segue o programa EES do roteiro~\cite{roteiro}, implementado em
Python, com as propriedades da água obtidas na biblioteca CoolProp: as do líquido e do vapor saturados em $T_{cond}$ e as da água de
resfriamento em $T_{resf}=(T_{a1}+T_{a2})/2$.

\subsection{Balanço de energia e região de condensação}

As vazões vêm dos volumes coletados, $\dot m=\rho V/\Delta t$, e as taxas de
calor são
\begin{equation}
  \dot Q_{resf}=\dot m_{resf}\,c_p\,(T_{a2}-T_{a1}),\qquad
  \dot Q_{cond}=\dot m_{cond}\,h_{lv}.
\end{equation}
A água de resfriamento recebe também o calor do subresfriamento do condensado.
Admitindo troca uniforme ao longo da serpentina, a condensação ocupa o
comprimento $L_{cond}=(\dot Q_{cond}/\dot Q_{resf})\,L_t$, com $L_t=1{,}43$ m.

\subsection{Coeficiente experimental}

Entre o vapor e a água de resfriamento há três resistências térmicas em série:
a do condensado ($R_i$), a da parede de vidro ($R_p$) e a da água ($R_e$):
\begin{equation}
  \dot Q_{cond}=\frac{T_{cond}-T_{resf}}{R_i+R_p+R_e},\quad
  R_i=\frac{1}{h_{cond}S_i},\quad
  R_p=\frac{\ln(d_{s,e}/d_{s,i})}{2\pi k_{v}L_{cond}},\quad
  R_e=\frac{1}{h_{resf}S_e},
  \label{eq:resist}
\end{equation}
com $S_i=\pi d_{s,i}L_{cond}$, $S_e=\pi d_{s,e}L_{cond}$ e
$k_v=\Kvidro$ W/(m\,K) (vidro Pyrex). Como o escoamento da água de resfriamento
é laminar ($\mathrm{Re}=\Reresf$), usou-se $h_{resf}=3{,}657\,k/d_{h,i}$
(parede a temperatura constante)~\cite{incropera}. Calculadas $R_p$ e $R_e$, a
Equação~\eqref{eq:resist} fornece $R_i$ e, daí, $h_{cond}$.

\subsection{Correlação de Dobson e Chato}

A correlação de Dobson e Chato~\cite{roteiro,incropera} é
\begin{equation}
  h_{cond}=0{,}023\,\mathrm{Re}_l^{0,8}\,\mathrm{Pr}_l^{0,4}
  \left(1+\frac{2{,}22}{\chi^{0,89}}\right)\frac{k_l}{d_{s,i}},\qquad
  \mathrm{Re}_l=\frac{4\dot m_{cond}(1-x)}{\pi d_{s,i}\mu_l},\quad
  \mathrm{Re}_v=\frac{4\dot m_{cond}\,x}{\pi d_{s,i}\mu_v},
  \label{eq:dc}
\end{equation}
em que $x$ é o título do vapor. Como $\mathrm{Re}_l$ e $\mathrm{Re}_v$ ficaram
abaixo de 2000, usou-se o parâmetro de Martinelli laminar,
$\chi_{ll}=[(1-x)/x]^{0,5}(\rho_v/\rho_l)^{0,5}(\mu_l/\mu_v)^{0,5}$.

O programa do roteiro avalia a correlação no título médio $x=0{,}5$. Como
melhoria, dividiu-se a região de condensação em $N=\N$ elementos, com o título
caindo linearmente de 1 (entrada) a 0 (fim da condensação): no elemento $i$,
$x_i=1-(i-\tfrac12)/N$. O coeficiente médio da região é a média dos
coeficientes locais, $\bar h=\frac1N\sum_{i=1}^{N}h(x_i)$, que equivale às
diretivas \textit{Duplicate} e \textit{Sum} do EES. Refinar a divisão até
$N=1000$ muda $\bar h$ em menos de 0,5\%.

As incertezas foram propagadas numericamente pelo método do GUM, a partir das
incertezas da tabela de medições, e são apresentadas como incertezas
expandidas (fator de abrangência 2).

\section{Resultados e discussão}

\begin{table}[H]
\centering
\caption{Medições}
\label{tab:medicoes}
\small
\begin{tabular}{lcc@{\hspace{2.2em}}lcc}
\toprule
\textbf{Grandeza} & \textbf{Valor} & \textbf{Incerteza} &
\textbf{Grandeza} & \textbf{Valor} & \textbf{Incerteza}\\
\midrule
$T_{cond}$ ($^{\circ}$C) & \Tcond & $\pm1\,^{\circ}$C & $V_{cond}$ (mL) & \Vcond & $\pm1\%$\\
$T_{a1}$ ($^{\circ}$C) & \Taum & $\pm1\,^{\circ}$C & $\Delta t_{cond}$ (s) & \tcond & $\pm5$ s\\
$T_{a2}$ ($^{\circ}$C) & \Tadois & $\pm1\,^{\circ}$C & $V_{AB}$ (V) & \VAB & $\pm5\%$\\
$V_{resf}$ (mL) & \Vresf & $\pm1\%$ & $I$ (A) & \Corrente & $\pm5\%$\\
$\Delta t_{resf}$ (s) & \tresf & $\pm5$ s & $T_{\mathrm{amb}}$ ($^{\circ}$C) & \Tamb & $\pm1\,^{\circ}$C\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[H]
\centering
\caption{Resultados}
\label{tab:resultados}
\small
\begin{tabular}{lcc}
\toprule
\textbf{Grandeza} & \textbf{Valor} & \textbf{Unidade}\\
\midrule
Vazão de água de resfriamento, $\dot m_{resf}$ & $\Mresf\pm\UMresf$ & g/s\\
Vazão de condensado, $\dot m_{cond}$ & $\Mcond\pm\UMcond$ & g/s\\
Potência elétrica, $V_{AB}I$ & $\Qel\pm\UQel$ & W\\
Calor retirado pela água de resfriamento, $\dot Q_{resf}$ & $\Qresf\pm\UQresf$ & W\\
Calor de condensação, $\dot Q_{cond}$ & $\Qcond\pm\UQcond$ & W\\
Comprimento da região de condensação, $L_{cond}$ & $\Lcond\pm\ULcond$ & m\\
\midrule
Coeficiente da água de resfriamento, $h_{resf}$ & $\Hresf\pm\UHresf$ & W/(m$^2$\,K)\\
Resistência da água, $R_e$ & $\Rext\pm\URext$ & K/W\\
Resistência da parede, $R_p$ & $\Rp\pm\URp$ & K/W\\
Resistência do condensado, $R_i$ & $\Ri\pm\URi$ & K/W\\
$h_{cond}$ experimental & \HexpNom$^{*}$ & W/(m$^2$\,K)\\
\midrule
$\mathrm{Re}_l$ e $\mathrm{Re}_v$ em $x=0{,}5$ & \Rel{} e \Rev & --\\
$h_{cond}$ Dobson e Chato, $x=0{,}5$ & $\Hdc\pm\UHdc$ & W/(m$^2$\,K)\\
$h_{cond}$ Dobson e Chato, média de $N=\N$ elementos & $\HdcN\pm\UHdcN$ & W/(m$^2$\,K)\\
\bottomrule
\end{tabular}

\vspace{0.3em}

\end{table}

\begin{figure}[H]
  \centering
  \includegraphics[width=\textwidth]{condensacao.png}
  \caption{(a) Resistências térmicas em série entre o vapor e a água de
  resfriamento, com a parcela de cada uma na resistência total e as barras de
  incerteza expandida. (b) Coeficiente de Dobson e Chato local nos $N=\N$
  elementos da região de condensação, sua média, o valor no título médio
  $x=0{,}5$ e o valor experimental.}
  \label{fig:condensacao}
\end{figure}

O vapor condensado corresponde a $\dot Q_{cond}=\Qcond$ W, cerca de
\FracEl\% da potência elétrica, mesma proporção observada na prática de
ebulição: o restante se perde pelo isolamento e pelas paredes do vaso. Como a
água de resfriamento recebeu $\Qresf$ W, a condensação ocupa cerca de
\Lcond{} m dos 1,43 m da serpentina, e o trecho final subresfria o condensado.

O valor experimental nominal, $h_{cond}=\HexpNom$ W/(m$^2$\,K), coincide com o
de Dobson e Chato no título médio, \HdcNom{} W/(m$^2$\,K). O cálculo local,
porém, mostra que o coeficiente não é uniforme
(Figura~\ref{fig:condensacao}b): ele cai nas duas pontas da região, perto de
$x=1$, onde quase não há líquido, e perto de $x=0$, onde o vapor já não
intensifica a troca. Por isso, o coeficiente no título médio fica acima da
média da região, que é \HdcNNom{} W/(m$^2$\,K), \DifN\% abaixo do valor em
$x=0{,}5$. Como o valor experimental representa toda a região de condensação,
é com essa média que ele deve ser comparado.

A comparação, no entanto, é apenas de ordem de grandeza. A resistência do
condensado, $R_i$, responde por só \FracRi\% da resistência total e é obtida
por diferença, enquanto a da água, $R_e$, responde por \FracRe\% e depende de
$L_{cond}$, isto é, de $\dot Q_{resf}$. Como $\dot Q_{resf}$ vem de uma
diferença de 9,5\,$^{\circ}$C medida com termopares de $\pm1\,^{\circ}$C e de
uma coleta de 60 s com $\pm5$ s, as medições da água de resfriamento respondem
por \FracVarResf\% da variância de $R_i$, cujo intervalo,
$(\Ri\pm\URi)$ K/W, inclui zero (Figura~\ref{fig:condensacao}a). Além disso,
$h_{resf}$ vem de um modelo laminar simplificado, e a correlação de Dobson e
Chato foi desenvolvida para escoamento anular turbulento, enquanto aqui
$\mathrm{Re}_l\approx\Rel$ e $\mathrm{Re}_v\approx\Rev$. Para reduzir a
incerteza, seria preciso medir melhor o aquecimento da água de resfriamento,
por exemplo com uma coleta mais longa e termômetros mais precisos.

Todo o código empregado nesta atividade está disposto no \href{}{Repositório GIT da Disciplina}

\begin{thebibliography}{9}
\bibitem{roteiro}
UFMG. \textit{EMA-103: Laboratório de Térmica -- Aulas 12 e 13: transferência
de calor por condensação}. Material didático da disciplina, 2026.

\bibitem{incropera}
BERGMAN, T. L.; LAVINE, A. S.; INCROPERA, F. P.; DEWITT, D. P.
\textit{Fundamentals of Heat and Mass Transfer}. 7. ed. Hoboken: Wiley, 2011.
\end{thebibliography}

\end{document}
