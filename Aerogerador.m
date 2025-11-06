clear all
close all
clc

%% PARAMETRIZAÇÃO %%
% Parametrização de aerogerador de 4.2MW

kv = 900; % [V/(rad/s)] Constante de velocidade
kt = kv; % [V/(rad/s)] Constante de torque
Ra = 100e-3; % [Ohm] Resistência de armadura 
La = 108e-3; % [H] Indutância de armadura
J = 9e6; % [kgm^2] Momento de inércia do rotor
R = 73.5; % [m] Raio do rotor
ro = 1.2; % [kg/m^3] Densidade do ar
b = 46e3; % [Nm/(rad/s)] Coeficiente de atrito total
A = pi*R^2; % [m^2] Área seção transversal do rotor
beta = 0; % [rad] Ângulo de pitch das pás

% Valores hipotéticos

Ea = 500;  % [V]
U = 10; % [m/s] 
w = 1;

%%
%---------------QUESTÃO 0.1--------------------------------------

% Intervalo de lambda

lambda = linspace(0, 15, 100);
betas = [0, 0.2, 0.5];
Cp = zeros(length(lambda), length(betas));


for i = 1:length(betas) % Calcular Cp para cada beta
    beta = betas(i);
    Cp(:, i) = (0.44 - 0.167*beta).*sin(pi*lambda./(15 - 0.3*beta)) - 0.16*lambda*beta;
end

figure;
plot(lambda, Cp, 'LineWidth', 1.8);
grid on;

xlabel('\lambda (razão de velocidade de ponta)');
ylabel('C_p (coeficiente de potência)');
title('Curvas de C_p(\lambda, \beta)');


legend(arrayfun(@(b) sprintf('\\beta = %.1f', b), betas, 'UniformOutput', false), ...
       'Location', 'best');


%%
%lambda = R*w/U;
%Cp = (0.44 - 0.167*beta)*sin((pi*lambda)/(15 - 0.3*beta)) - 0.16*lambda*beta;
funCp = @(lambda, beta)((0.44 - 0.167*beta)*sin((pi*lambda)/(15 - 0.3*beta)) - 0.16*lambda*beta);

Pw = 0.5*ro*A*U^3;

% Pm = tau*w = Cp*Pw
% tau = Cp*Pw/w
% tau = (Cp * Pw)/w;

funTau = @(w) ((funCp(R*w/U, beta) * Pw)/w);

% ./w  ao inves de /w faria multiplicação elemento a elemento? perguntar
% pro lelis
% .* é o equivalente para multiplicação


%%  Equilíbrio %% 
% Feito substituindo uma edo na outra e removendo dependencia de ia em w
w0 = 1;   % CI 
% se CI = 0 da erro porque tau/0

w = fzero(@(w) funTeste(w, funTau, b, kv, Ra, Ea), w0);

ia = (kv*w - Ea)/Ra;






% DINAMICA
% [tvec,yvec] = ode45(@(t, y) funOmega(J, tau, b, y, kv, ia),[0 10*60],[0 0]);
% 
% 
% [t2vec,xvec] = ode45(@(t, x) funIa(La, kv, w, Ra, Ea, x),[0 10*60],[0 0]);
% 
% iavec = xvec(:, 1);

%% Item 4 %%

%Ppe = Ra*ia^2; % Perdas elétricas

%Ppm = b*w^2; % Perdas mecânicas

%Pg = Ea*ia; % Potência gerada

%Pe = tau*w; % potencia injetada?




