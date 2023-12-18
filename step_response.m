clear all; close all; clc;

%% Inserimento valori
m = 2500;       %% massa
b = 350;        %% smorzamento
k = 5*10^5;     %% cost. elastica lineare
ze = 0.35;      %% stato iniziale z

figure;
hold on; box on; zoom on; grid on;

s = tf('s');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% definire FdT
G = 1/(m*s^2+s*b+k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% risposta impulsiva fdt
time = 0:0.01:100;

YY_i = impulse(G, time);
YY_s = step(G, time);

% plot
plot(time, YY_i, 'r', 'DisplayName', 'impulse response G(s)', 'LineWidth', 1.3);
plot(time, YY_s, 'b', 'DisplayName', 'step response G(s)', 'LineWidth', 1.3);

title('Risposta fdt')
xlabel('tempo [s]')
ylabel('posizione')
legend(["impulse"; "step"]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

