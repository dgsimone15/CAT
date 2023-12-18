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
%%% fdt
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


%% uscita fdt ad ingresso uu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modello funzione di trasferimento
[YY, TT] = lsim(G,uu, tt);

%stampa della fuznione di trasferimento
figure(2);
hold on; grid on; zoom on;
plot(TT,YY);
%%plot(TT, uu); //stampa ingresso
title('Traiettoria di stato y(t) - calcolata dalla fdt')
xlim([0, 10])
xlabel('tempo [s]')
ylabel('posizione')
legend('y(t) - posizione')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% risposta impulsiva fdt
time = 0:0.01:100;

figure(3);
hold on; grid on; zoom on;
YY_i = impulse(G, time);
YY_s = step(G, time);

% plot
plot(time, YY_i, 'r', 'DisplayName', 'impulse response G(s)', 'LineWidth', 1.3);
plot(time, YY_s, 'b', 'DisplayName', 'step response G(s)', 'LineWidth', 1.3);

title('Risposta fdt')
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

%% controllo fdt
Gcontrol = tf(modello) %%ritorna la funzione di trasferimento dal modello, utile per confronto con la tf calcolata
Id = eye(2);

%% stampa confronto uscite tra fdt e modello LTI
figure;
hold on; box on; zoom on; grid on;
plot(TT,YY, 'r') %%fdt
plot(TT,YY2, 'g') %%LTI
title('Confronto modello LTI e fdt')
legend('y(t) fdt', 'y(t) - LTI');
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