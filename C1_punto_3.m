clear all; close all; clc

%solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-3;
omega_plot_max = 1e6;

%% Parametri

m = 2500;       %% mass
b = 350;        %% smorzamento
k = 5*10^5;     %% cost. elastica lineare
ze = 0.35;      %% stato iniziale z

beta = 0.1;     %% coeff. di non lienarità elastica
n = 3;          %% coeff. di non lienarità elastica

%% Calcolo coppia di equilibrio
x_1e = ze;

x_2e = 0;
u_e  = 9.81*m+k*x_1e;
x_e  = [x_1e;x_2e];


%% Linearizzazione e calcolo FdT
% matrici del sistema
A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;

s  = tf('s');
GG = C*inv(s*eye(2) - A)*B + D

%% Specifiche

% errore a regime
WW = 1;
DD = 1;
e_star = 0.02;

% attenuazione disturbo sull'uscita
A_d = 30;
omega_d_min = omega_plot_min;
omega_d_MAX = 0.05;

% attenuazione disturbo di misura
A_n = 75;
omega_n_min = 8e4;
omega_n_MAX = 8e6;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_star = 5;
T_star = 0.01;

% Margine di fase
Mf_esp = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% progettare regolatore statico RR_s

mu_s_error = (WW + DD)/e_star
% calcolo del guadagno voluto per avere errore e_star

mu_s_RR_s = mu_s_error/abs(evalfr(GG,j*0));
% calcolo guadagno del regolatore statico in base al guadagno voluto/guadagno di GG a 0

mu_s_d = 10^(A_d/20)/abs(evalfr(GG,j*omega_d_MAX))
% calcolo guadagno del regolatore statico all'infinito (alla pulsazione massima per errore d)

RR_s = max(mu_s_RR_s, mu_s_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sistema esteso
GG_e = RR_s*GG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mappare specifiche sul diagramma di Bode

figure(1);
hold on;

%%%%%%%%%%%
%%% mappare attenuazione disturbo uscita (in diagramma di ampiezza)

% coordinate orizzontali
patch_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];

% coordinate verticali
patch_d_y = [A_d; A_d; -200; -200];

patch(patch_d_x, patch_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
%%%%%%%%%%%

%%%%%%%%%%%
%%% mappare attenuazione disturbo di misura (in diagramma di ampiezza)

% coordinate orizzontali
patch_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];

% coordinate verticali
patch_n_y = [-A_n; -A_n; 200; 200];

patch(patch_n_x, patch_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
%%%%%%%%%%%

%%%%%%%%%%%
%%% mappare tempo di assestamento (in diagramma di ampiezza)

% calcolo di xi_star a partire dalla sovraelongazione percentuale
xi_star = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);

% calcolo margine di fase
M_f_S = 100*xi_star;

%trovare il margine di fase maggiore per soddisfare sia quello della specifica sia quello trovato dalla sovraelongazione
M_f = max(M_f_S,Mf_esp);

% pulsazione minima per garanitre tempo di assestamento
omega_c_min = 460/(M_f*T_star);

% coordinate orizzontali
patch_omega_c_x = [omega_plot_min; omega_c_min; omega_c_min; omega_plot_min];

% coordinate verticali
patch_omega_c_y = [0; 0; -200; -200];

patch(patch_omega_c_x, patch_omega_c_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
%%%%%%%%%%%

%Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;

Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

%%%%%%%%%%%
%%% mappare margine di fase (in diagramma di fase)

% creiamo pulsazione massima
omega_c_max = omega_n_min;

% coordinate orizzontali
patch_M_f_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];

% coordinate verticali
patch_M_f_y = [M_f-180; M_f-180; -270; -270];

patch(patch_M_f_x, patch_M_f_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

Legend_arg = ["G(j\omega)"; "M_f"]; % Legenda colori
legend(Legend_arg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;