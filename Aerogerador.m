clear all
close all
clc

%% PARAMETRIZA√á√ÉO %%
% Parametriza√ß√£o de aerogerador de 4.2MW

kv = 900; % [V/(rad/s)] Constante de velocidade
kt = kv; % [V/(rad/s)] Constante de torque
Ra = 100e-3; % [Ohm] Resist√™ncia de armadura 
La = 108e-3; % [H] Indut√¢ncia de armadura
J = 9e6; % [kgm^2] Momento de in√©rcia do rotor
R = 73.5; % [m] Raio do rotor
rho = 1.2; % [kg/m^3] Densidade do ar
b = 46e3; % [Nm/(rad/s)] Coeficiente de atrito total
A = pi*R^2; % [m^2] √?rea se√ß√£o transversal do rotor
beta = 0; % [rad] √Çngulo de pitch das p√°s

% Valores hipot√©ticos

Ea = 500;  % [V]
U = 10; % [m/s] 
%%
%---------------QUEST√O 0.1--------------------------------------

% Intervalo de lambda

lambda = linspace(0, 15, 100);
betas = [0, 0.2, 0.5];
Cp = zeros(length(lambda), length(betas));


for i = 1:length(betas) % Calcular Cp para cada beta
    beta = betas(i);
    Cp(:, i) = (0.44 - 0.167*beta).*sin(pi*lambda./(15 - 0.3*beta)) - 0.16*lambda*beta;
end

figure;
plot(lambda, Cp, 'LineWidth', 1);
grid on;

xlabel('\lambda (raz„o de velocidade de ponta)');
ylabel('C_p (coeficiente de potÍncia)');
title('Curvas de C_p(\lambda, \beta)');


legend(arrayfun(@(b) sprintf('\\beta = %.1f', b), betas, 'UniformOutput', false), ...
       'Location', 'best');


%%
%lambda = R*w/U;
%Cp = (0.44 - 0.167*beta)*sin((pi*lambda)/(15 - 0.3*beta)) - 0.16*lambda*beta;
funCp = @(lambda, beta)((0.44 - 0.167*beta)*sin((pi*lambda)/(15 - 0.3*beta)) - 0.16*lambda*beta);

Pw = 0.5*rho*A*U^3;

% Pm = tau*w = Cp*Pw
% tau = Cp*Pw/w
% tau = (Cp * Pw)/w;

funTau = @(w) ((funCp(R*w/U, beta) * Pw)/w);

% ./w  ao inves de /w faria multiplica√ß√£o elemento a elemento? perguntar
% pro lelis
% .* √© o equivalente para multiplica√ß√£o


%%  Equil√≠brio %% 
% Feito substituindo uma edo na outra e removendo dependencia de ia em w
w0 = 1;   % CI 
% se CI = 0 da erro porque tau/0

w = fzero(@(w) funTeste(w, funTau, b, kv, Ra, Ea), w0);

ia = (kv*w - Ea)/Ra;





%% 
%-----------------------QUEST√O 4----------------------------------

Ppe = Ra*ia^2; % Perdas el√©tricas

Ppm = b*w^2; % Perdas mec√¢nicas

Pg = Ea*ia; % Pot√™ncia gerada

Pe = funTau(w)*w; % PotÍncia injetada

%% 
%-----------------------QUEST√O 5----------------------------------

x0 = [10; 2; 400; 10000];
x = fsolve(@(x) funNominal(x, rho, A, b, kv, R, Ra), x0);

U_nom = x(1, :);
w_nom = x(2, :);
Ea_nom = x(3, :);
ia_nom = x(4, :);


%%
%-----------------------QUEST√O 9----------------------------------
% DINAMICA

[tvec, xvec] = ode45(@(t,x) funDin_aerogerador(J, funTau, b, La, Ra, Ea, kv, x), [0 10], [0 1]);
iavec = xvec(:,1);
wvec = xvec(:,2);

% Plot din‚mica

subplot(2,1,1)
plot(tvec, iavec, 'r', 'LineWidth', 1)
xlabel('Tempo [s]')
ylabel('Corrente [A]')
grid on

subplot(2,1,2)
plot(tvec, wvec, 'b', 'LineWidth', 1)
xlabel('Tempo [s]')
ylabel('Velocidade de rotaÁ„o [rad/s]')
grid on

sgtitle('Din‚mica das Vari·veis de Estado', 'FontWeight', 'bold', ...
    'FontSize', 12)






