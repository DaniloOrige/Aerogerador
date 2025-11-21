clear all
close all
clc

%% PARAMETRIZA√á√ÉO %%
% Parametriza√ß√£o de aerogerador de 4.2MW

kv = 900*2; % [V/(rad/s)] Constante de velocidade
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

beta0 = 0; %[graus]
beta1 = 15; %[graus]
beta2 = 30; %[graus]

% PESQUISAR BETAS QUE S√O USADOS NORMALMENTE

lambda = linspace(0, 15, 100);
betas = [0, beta1*pi/180, beta2*pi/180];
Cp = zeros(length(lambda), length(betas));


for i = 1:length(betas) % Calcular Cp para cada beta
    beta = betas(i);
    Cp(:, i) = (0.44 - 0.167*beta).*sin(pi*lambda./(15 - 0.3*beta)) - 0.16*lambda*beta;
end

Cpmax = max(Cp,[], 1); % Extraindo valor m·ximo para cada beta

Cpmax0 = Cpmax(:, 1);
Cpmax15 = Cpmax(:, 2);
Cpmax30 = Cpmax(:, 3);


figure;
plot(lambda, Cp, 'LineWidth', 1.5);
grid on;

xlabel('\lambda (raz„o de velocidade de ponta)');
ylabel('C_p (coeficiente de potÍncia)');
title('Curvas de C_p(\lambda, \beta)');


legend(arrayfun(@(b) sprintf('\\beta = %.1f', b), betas*180/pi, 'UniformOutput', false), ...
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

%%
%-----------------------QUEST√O 2----------------------------------
% Feito substituindo uma edo na outra e removendo dependencia de ia em w
w0 = 1;   % CI 
% se CI = 0 da erro porque tau/0

w = fzero(@(w) funTeste(w, funTau, b, kv, Ra, Ea), w0);

ia = (kv*w - Ea)/Ra;


%% 
%-----------------------QUEST√O 4----------------------------------

Ppe = Ra*ia^2; % Perdas el√©tricas

Ppm = b*w^2; % Perdas mec√¢nicas

Pg = kv*ia*w; % Pot√™ncia gerada

Pe = Pg - Ppe - Ppm; % PotÍncia entregue

%% 
%-----------------------QUEST√O 5----------------------------------

lambdaopt = 7.5;
Pnom = 4.2*10^6;
x0 = [60; 5000];

options = optimset('MaxFunEvals',1000,'Display','none');
[x,fval] = fsolve(@(x) funNominal(x, rho, A, b, kv, R, Ra, Pnom, lambdaopt, Cpmax0), x0,options);

U_nom = x(1, :);
ia_nom  = x(2, :);
Ea_nom = Pnom/ia_nom;
w_nom = (U_nom*lambdaopt)/R;


%%
%-----------------------QUEST√O 6----------------------------------

Uvec_menor = 3:U_nom;

Ea_vec = zeros(size(Uvec_menor));
wvec_menor  = zeros(size(Uvec_menor));

options = optimset('MaxFunEvals',1000,'Display','none');
x0 = [0; 0];

for k = 1:length(Uvec_menor)
    
    U = Uvec_menor(k);

    [x, fval] = fsolve(@(x) funMenor(x, rho, A, b, kv, R, Ra, U, lambdaopt, Cpmax0), x0, options);

    Ea_vec(k) = x(1);
    wvec_menor(k)  = (U*lambdaopt)/R;

  
end

% Plot 1
figure;

subplot(2,1,1);
plot(Uvec_menor, Ea_vec, 'g', 'LineWidth', 1.5);
xlabel('Velocidade do vento U [m/s]');
ylabel('Tens„o de armadura E_a [V]');
grid on;

subplot(2,1,2);
plot(Uvec_menor, wvec_menor, 'r', 'LineWidth', 1.5);
xlabel('Velocidade do vento U [m/s]');
ylabel('Velocidade angular \omega [rad/s]');
grid on;

sgtitle('Comportamento do Aerogerador na Regi„o Subnominal (U < Unom)', ...
    'FontWeight','bold', 'FontSize',12);

%%
% U > U_nom

Uvec_maior = U_nom:20;
Beta_vec = zeros(size(Uvec_maior));
wvec_maior  = zeros(size(Uvec_maior));

for k = 1:length(Uvec_maior)

    U = Uvec_maior(k);
    Beta = fzero(@(beta) funMaior(lambdaopt, U, R, rho, A, ia_nom, kv, beta, b), 0);
    Cp_teste = (0.44 - 0.167*Beta)*sin((pi*lambdaopt)/(15 - 0.3*Beta)) - 0.16*lambdaopt*Beta;

    Beta_vec(k) = Beta;
    wvec_maior(k)  = (U*lambdaopt)/R;

end

% Plot 2

figure;
subplot(2,1,1);
plot(Uvec_maior, Beta_vec*180/pi, 'LineWidth', 1.5)
xlabel('Velocidade do vento U [m/s]');
ylabel('\beta  [graus]');
grid on;

subplot(2,1,2);
plot(Uvec_maior, wvec_maior, 'r', 'LineWidth', 1.5);
xlabel('Velocidade do vento U [m/s]');
ylabel('Velocidade angular \omega [rad/s]');
grid on;

sgtitle('Comportamento do Aerogerador na Regi„o Subnominal (U > Unom)', ...
    'FontWeight','bold', 'FontSize',12);

%%

%-----------------------QUEST√O 9----------------------------------
% DINAMICA

[tvec, xvec] = ode45(@(t,x) funDin_aerogerador(J, funTau, b, La, Ra, Ea, kv, x), [0 10], [0 1]);
iavec = xvec(:,1);
wvec = xvec(:,2);

% Plot din‚mica

subplot(2,1,1)
plot(tvec, iavec, 'r', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Corrente [A]')
grid on

subplot(2,1,2)
plot(tvec, wvec, 'b', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade de rotaÁ„o [rad/s]')
grid on

sgtitle('Din‚mica das Vari·veis de Estado', 'FontWeight', 'bold', ...
    'FontSize', 12)

%%
%-----------------------QUEST√O 10----------------------------------


%%












