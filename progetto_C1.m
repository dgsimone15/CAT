close all; clc;
clear all; %calo significativo delle performance, attivare solo se serve

%% parametri sistema
s = tf('s'); %%variabile della funzione di trasferimento
tt = 0:0.01:10; % tempo, da impostare ad almeno -> 0 a 10 secondi con passo 0.01, per risultati affidabili

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inserimento valori
m = 2500;       %% massa
b = 350;        %% smorzamento
k = 5*10^5;     %% cost. elastica lineare
ze = 0.35;      %% stato iniziale z

beta = 0.1;     %% coeff. di non lienarità elastica
n = 3;          %% coeff. di non lienarità elastica

% condizioni iniziali
pos_init = ze; % [m]
vel_init = 0; % [m/s]
x0 = [pos_init; vel_init];

%%uu(1:length(tt),1) = 500; ingresso di test con valore costante
uu(1:length(tt),1) = 9.81*m+k*ze; %%u(t) del sistema, ingresso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parte 1

% matrici del sistema linearizzato
A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FDT
G = C/(s*eye(2) - A)*B + D

%%% plot diagramma di Bode
figure(1);
hold on;
legend(["G(j\omega)"])
bode(G);
title("FDT - G")
legend(["G(j\omega)"])
grid on; zoom on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% uscita FDT ad ingresso uu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modello funzione di trasferimento
[YY, TT] = lsim(G,uu, tt);

%stampa della fuznione di trasferimento
figure(2);
hold on; grid on; zoom on;
plot(TT,YY);
%%plot(TT, uu); //stampa ingresso
title('Traiettoria di stato y(t) - calcolata dalla FDT')
xlim([0, 10])
xlabel('tempo [s]')
ylabel('posizione')
legend('y(t) - posizione')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% risposta impulsiva FDT
time = 0:0.01:100;

figure(3);
hold on; grid on; zoom on;
YY_i = impulse(G, time);
YY_s = step(G, time);

% plot
plot(time, YY_i, 'r', 'DisplayName', 'impulse response G(s)', 'LineWidth', 1.3);
plot(time, YY_s, 'b', 'DisplayName', 'step response G(s)', 'LineWidth', 1.3);

title('Risposta FDT')
xlabel('tempo [s]')
ylabel('posizione')
legend(["impulse"; "step"]);

%% test
%{

%% G fattorizzato
G_manual = 1/(m*s^2+s*b+k);
G_fatt = zpk(G) %%zero pole gain model -> forma fattorizzata del primo tipo, p->1/m = 0.0004


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% modello LTI

% state-space model e soluzione dell'equazione
modello = ss(A, B, C, D);

[YY2, TT, XX] = lsim(modello, uu, tt, x0); %%simulazione modello dinamico

% stampa stato del sistema LTI
figure;
hold on; box on; zoom on; grid on;
plot(TT, XX);
%%plot(TT, uu); //stampa ingresso
title('Traiettoria di stato y(t) - Modello LTI')
xlim([0, 10]);
xlabel('tempo [s]')
ylabel('stato')
legend('y(t) - posizione', 'x2 - velocità');

% Dai risultati ottenuti è possibile osservare come nel momento in cui
% l'oggetto comincia a risalire, la velocità passa da valore negativo a
% valore positivo.

%% controllo FDT
Gcontrol = tf(modello) %%ritorna la funzione di trasferimento dal modello, utile per confronto con la tf calcolata
Id = eye(2);

%% stampa confronto uscite tra FDT e modello LTI
figure;
hold on; box on; zoom on; grid on;
plot(TT,YY, 'r') %%FDT
plot(TT,YY2, 'g') %%LTI
title('Confronto modello LTI e FDT')
legend('y(t) FDT', 'y(t) - LTI');
%%note:
%sovrapposizione perfetta per: x0 = (0,0)
% con x0 != 0 la funzione LTI ha un comportamento diverso perché considera
% l'origine del sistema

%% stampa margine di fase funzione di trasferimento
figure;
margin(G); % privo di senso dal momento che la funzione parte da valori nettamente inferiori allo 0db giá nell'origine

%}


%% animazioni

%{
%%
figure
tt = 0:1e-3:0.3; % intervallo temporale: da 0 ad 0.3 con passo 0.001




Y = [0, 5, 0:((+9 +55-4)/9):70-4, 55-4, 55-2]
length(Y)

    plot([0 60],[0 0],'k','LineWidth',4) % ground.
        hold on; box on; zoom on; grid on;
    axis([0 50  0 150])    
    plot(20+[0 0 0 1.5 -1.5 1.5 -1.5 1.5 -1.5 1.5 -1.5 0 0 0],[0, 5, 5:((+9 +55-4)/9):70-4, 70-4, 70-2],'r','LineWidth',2) % spring

%%piatto
    patch ([10; 40; 40; 10],[70-2; 70-2; 70-2+15; 70-2+15] ,'r','FaceAlpha',0.1) ; % rettangolo rosso con vertici definiti da x ed y
%%smorzatore
    patch ([27.5; 32.5; 32.5; 27.5],[30; 30; 40; 40] ,'v','FaceAlpha',0.1) ; % rettangolo rosso con vertici definiti da x ed y
    plot([30 30],[0 30],'k','LineWidth',4) % ground
    plot([30 30],[40 70-2],'k','LineWidth',4) % ground
%%u
    patch ([24; 26; 26; 24],[10; 10; 50; 50] ,'v','FaceAlpha',0.1) ; % rettangolo rosso con vertici definiti da x ed y
    plot([25 25],[0 10],'k','LineWidth',4) % ground.
    plot([25 25],[50 70-2],'k','LineWidth',4) % ground


    hold on; box on; zoom on; grid on;
    uu = 9.81*m + k*[ze;0];
    yy = lsim(G, uu, tt);

%}

%%%%%%%%%%%%%%%%%%%%%%% PUNTO 3 %%%%%%%%%%%%%%%%%%%%%%%

%solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-3;
omega_plot_max = 1e6;

%% Calcolo coppia di equilibrio
x_1e = ze;
x_2e = 0;
x_e  = [x_1e;x_2e];

u_e  = 9.81*m+k*x_1e;



%% Linearizzazione e calcolo FdT
% matrici del sistema
A = [0, 1; -k/m, -b/m];
B = [0; 1/m];
C = [1, 0];
D = 0;

GG = C/(s*eye(2) - A)*B + D

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

% Sovraelongazione massima e tempo d'assestamento al 5%
S_star = 5;
epsilon = 5;
T_star = 0.01;

% Margine di fase
Mf_esp = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% progettare REGOLATORE STATICO RR_s

mu_s_error = (WW + DD)/e_star
% calcolo del guadagno voluto per avere errore e_star

mu_s_RR_s = mu_s_error/abs(evalfr(GG,j*0));
% calcolo guadagno del regolatore statico come rapporto tra guadagno voluto e guadagno di GG a 0

mu_s_d = 10^(A_d/20)/abs(evalfr(GG,j*omega_d_MAX));
% calcolo guadagno del regolatore statico all'infinito (alla pulsazione massima per errore d)

RR_s = max(mu_s_RR_s, mu_s_d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sistema esteso
GG_e = RR_s*GG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mappare specifiche sul diagramma di Bode

figure(4);
hold on;

%%%%%%%%%%%
%%% mappare attenuazione disturbo uscita (in diagramma di ampiezza)

% coordinate orizzontali
patch_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];

% coordinate verticali
patch_d_y = [A_d; A_d; -200; -200];

patch(patch_d_x, patch_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
%%%%%%%%%%%
%%
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

Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

%Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Design del REGOLATORE DINAMICO RR_d

% ci troviamo in uno scenario B che risolviamo tramite una rete anticipatrice

%determinare una pulsazione che si vuole prendere nella zona in cui si rispetta il vincolo di M_f
omega_c_star = 900;

%calcoliamo M_star e phi_star
M_omega_c_star = abs(evalfr(GG_e,j*omega_c_star));
arg_omega_c_star = rad2deg(angle(evalfr(GG_e,j*omega_c_star)));

M_star = 1/M_omega_c_star;
phi_star = M_f + 5 - 180 - arg_omega_c_star; % prendendo un piccolo margine in più di +5

% FORMULE DI INVERSIONE (calcoliamo tau e alpha_tau)
tau = (M_star - cosd(phi_star))/(omega_c_star*sind(phi_star));
alpha_tau = (cosd(phi_star) - 1/M_star)/(omega_c_star*sind(phi_star));

if(min(tau, alpha_tau) < 0)
    fprintf('Valori negativi. non corretto!\n');
    return;
end

RR_d = (1+tau*s)/(1+alpha_tau*s);

% regolatore in cui inseriamo un polo in alta frequenza (prima di omega_n_min)
% in questo caso il polo è stato scelto a 15000 (<80000)
R_high_frequency = 1/(1 + s/15000);

% il polo inserito influenza leggermente la fase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = RR_d*RR_s*R_high_frequency;
%RR = RR_d*RR_s; %funzione di controllo priva di polo in alta frequenza

LL = RR*GG; % funzione di anello aperto


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% verificare che la funzione d'anello LL non violi le zone proibite
figure(5);
hold on;

patch(patch_d_x, patch_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

patch(patch_n_x, patch_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

patch(patch_omega_c_x, patch_omega_c_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
%%TODO correggere legenda
legend(Legend_mag);
margin(LL,{omega_plot_min,omega_plot_max}); % Plot Bode con margini di stabilità
grid on; zoom on;

% margine di fase
patch(patch_M_f_x, patch_M_f_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

legend(Legend_arg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mostrare che la risposta al gradino rispetti le specifiche

T_simulation = 2*T_star; %(+1000 per vedere errore a regime)
FF = LL/(1+LL);

%%TEST FF_s -> risposta dell'anello chiuso privo di regolatore dinamico
%{
    FF_s = GG_e/(1+GG_e);
    figure();
    hold on;
    step(WW*FF_s, T_simulation);
%}


%risposta a gradino di tutto il sistema retroazionato
[y_step, t_step] = step(WW*FF, T_simulation);

figure(6);
hold on;
plot(t_step, y_step);

%
LV = evalfr(WW*FF, 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sovraelongazione
patch_s_x = [0; T_simulation; T_simulation; 0];
patch_s_y = [LV*(1+S_star/100); LV*(1+S_star/100); 2*LV; 2*LV];

patch(patch_s_x,patch_s_y,'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tempo di assestamento
patch_T_x = [T_star; T_simulation; T_simulation; T_star];
patch_T_y = [LV*(1+epsilon/100); LV*(1+epsilon/100); 2*LV; 2*LV];
patch(patch_T_x,patch_T_y,'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vieta valori più piccoli dell'1% rispetto al valore asintotico
patch_T_y_down = [LV*(1-epsilon/100); LV*(1-epsilon/100); 0; 0];

patch(patch_T_x,patch_T_y_down,'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PUNTO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test sul sistema linearizzato con w, d e n dati da testo

WW = 1;
DD = 0;
NN = 0;

SS = 1/(1+LL); % funzione di sensitività
% FF -> unzione di sensitività complementare

figure(7)
bode(SS);
grid on; zoom on;
title("FDT SS");

figure(8);
bode(FF);
grid on; zoom on;
title("FDT FF")

%{
usiamo un unico tt per permettere di stampare la YY generale, é un suicidio
prestazionale!!
%}
%% Simulazione disturbo in uscita
tt = 0:1e-3:2e3;
for k=1:1:4
    DD = DD + sin(0.01*k*tt);
end
y_d = lsim(SS,DD,tt);
% plot
figure(10)
hold on; grid on; zoom on;
title("y_d - disturbo d'uscita");
plot(tt,y_d,'r');
plot(tt,DD,'g');


%% Simulazione disturbo di misura
for k=1:1:4
    NN = NN + 0.2*sin((8*1e4)*k*tt);
end
y_n = lsim(FF,NN,tt);
% plot
figure(9)
hold on; grid on; zoom on;
title("y_n - disturbo lettura");
xlim([0; 1]);
plot(tt,y_n,'r');
plot(tt,NN,'g');


%% risposta al gradino del sistema globale
[y_step, t_step] = step(WW*FF, tt);
YY = y_step+y_d+y_n;
% plot
figure(11)
plot(tt,YY);
xlim([0;1]);
title("YY - risposta al gradino dell'intero sistema")