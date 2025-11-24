clear all
close all
clc

%% PARAMETRIZAÃ‡ÃƒO %%
% ParametrizaÃ§Ã£o de aerogerador de 4.2MW

kv = 900*2; % [V/(rad/s)] Constante de velocidade
kt = kv; % [V/(rad/s)] Constante de torque
Ra = 100e-3; % [Ohm] ResistÃªncia de armadura 
La = 108e-3; % [H] IndutÃ¢ncia de armadura
J = 9e6; % [kgm^2] Momento de inÃ©rcia do rotor
R = 73.5; % [m] Raio do rotor
rho = 1.2; % [kg/m^3] Densidade do ar
b = 46e3; % [Nm/(rad/s)] Coeficiente de atrito total
A = pi*R^2; % [m^2] Ã?rea seÃ§Ã£o transversal do rotor
beta = 0; % [rad] Ã‚ngulo de pitch das pÃ¡s

% Valores hipotÃ©ticos

Ea = 500;  % [V]
U = 10; % [m/s] 
%%
%---------------QUESTÃO 0.1--------------------------------------

% Intervalo de lambda

beta0 = 0; %[graus]
beta1 = 15; %[graus]
beta2 = 30; %[graus]

% PESQUISAR BETAS QUE SÃO USADOS NORMALMENTE

lambda = linspace(0, 15, 100);
betas = [0, beta1*pi/180, beta2*pi/180];
Cp = zeros(length(lambda), length(betas));


for i = 1:length(betas) % Calcular Cp para cada beta
    beta = betas(i);
    Cp(:, i) = (0.44 - 0.167*beta).*sin(pi*lambda./(15 - 0.3*beta)) - 0.16*lambda*beta;
end

Cpmax = max(Cp,[], 1); % Extraindo valor máximo para cada beta

Cpmax0 = Cpmax(:, 1);
Cpmax15 = Cpmax(:, 2);
Cpmax30 = Cpmax(:, 3);


figure;
plot(lambda, Cp, 'LineWidth', 1.5);
grid on;

xlabel('\lambda (razão de velocidade de ponta)');
ylabel('C_p (coeficiente de potência)');
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

% ./w  ao inves de /w faria multiplicaÃ§Ã£o elemento a elemento? perguntar
% pro lelis
% .* Ã© o equivalente para multiplicaÃ§Ã£o

%%
%-----------------------QUESTÃO 2----------------------------------
% Feito substituindo uma edo na outra e removendo dependencia de ia em w
w0 = 1;   % CI 
% se CI = 0 da erro porque tau/0

w = fzero(@(w) funTeste(w, funTau, b, kv, Ra, Ea), w0);

ia = (kv*w - Ea)/Ra;


%% 
%-----------------------QUESTÃO 4----------------------------------

Ppe = Ra*ia^2; % Perdas elÃ©tricas

Ppm = b*w^2; % Perdas mecÃ¢nicas

Pg = kv*ia*w; % PotÃªncia gerada

Pe = Pg - Ppe - Ppm; % Potência entregue

%% 
%-----------------------QUESTÃO 5----------------------------------

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
%-----------------------QUESTÃO 6----------------------------------

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
ylabel('Tensão de armadura E_a [V]');
grid on;

subplot(2,1,2);
plot(Uvec_menor, wvec_menor, 'r', 'LineWidth', 1.5);
xlabel('Velocidade do vento U [m/s]');
ylabel('Velocidade angular \omega [rad/s]');
grid on;

sgtitle('Comportamento do Aerogerador na Região Subnominal (U < Unom)', ...
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

sgtitle('Comportamento do Aerogerador na Região Subnominal (U > Unom)', ...
    'FontWeight','bold', 'FontSize',12);

%%

%-----------------------QUESTÃO 9----------------------------------
% DINAMICA

[tvec, xvec] = ode45(@(t,x) funDin_aerogerador(J, funTau, b, La, Ra, Ea, kv, x), [0 10], [0 1]);
iavec = xvec(:,1);
wvec = xvec(:,2);

% Plot dinâmica

subplot(2,1,1)
plot(tvec, iavec, 'r', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Corrente [A]')
grid on

subplot(2,1,2)
plot(tvec, wvec, 'b', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade de rotação [rad/s]')
grid on

sgtitle('Dinâmica das Variáveis de Estado', 'FontWeight', 'bold', ...
    'FontSize', 12)

%%
%-----------------------QUESTÃO 10----------------------------------
% 

w_10a = out.questao10a;
w_ref10a = out.wref_10a;
u_10a = out.controle10a;
q_10a= out.perturbacao10a;


w_10b = out.questao10b;
w_ref10b = out.w_ref10b;
u_10b = out.controle10b;
q_10b = out.perturbacao10b;

t_simulink = out.tout;


% Plot questão 10a


figure;
subplot(3,1,1)
plot(t_simulink, w_10a, 'g', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade angular das pás [rad/s]')
grid on
hold on
plot(t_simulink, w_ref10a, 'Color' ,'#FF0000', 'LineWidth', 2.0, 'LineStyle','--');
legend('Resposta do sistema', 'Refência')

subplot(3,1,2)
plot(t_simulink, u_10a, 'b', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Entrada de controle')
grid on

subplot(3,1,3)
plot(t_simulink, q_10a, 'm', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade do vento [m/s]')
grid on
sgtitle('Ensaio resposta ao degrau com U = Unominal e referência variável', 'FontWeight', 'bold', ...
    'FontSize', 12)


figure;
subplot(3,1,1)
plot(t_simulink, w_10b, 'Color', '#F1C338', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade angular das pás [rad/s]')
grid on
hold on
plot(t_simulink, w_ref10b,'Color' ,'#FF00FF', 'LineWidth', 2.0, 'LineStyle','--');
legend('Resposta do sistema', 'Refência')

subplot(3,1,2)
plot(t_simulink, u_10b, 'b', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Entrada de controle')
grid on

subplot(3,1,3)
plot(t_simulink, q_10b, 'm', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade do vento [m/s]')
grid on



sgtitle('Ensaio de resposta ao degrau com referência', 'FontWeight', 'bold', ...
    'FontSize', 12)
%%
%-----------------------QUESTÃO 11----------------------------------
% 
% Uvec = (3:20)';
% vectempo = (0:4.5:76.5)';
% 
% dadosts = timeseries(vectempo, Uvec);
% dados_ts.UserData.Interpolation = 'zoh';
% 
% stairs(vectempo, Uvec, 'r', 'LineWidth', 1.5)
% 
% xlabel('Tempo [s]')
% ylabel('Vento [m/s]')
% grid on

vecU11 = linspace(3, 20, 18);
vecU11 = [3, vecU11];
vecstepinterval = zeros(1, length(vecU11));

for k = 1:length(vecU11) - 1
    vecstepinterval(k + 1) = 8*k;
end

w11 = out.w11;
tout = out.tout;
Ea11 = out.Ea11;
beta11 = out.beta11;
q11 = out.perturbacao11;
w_ref11 = out.w_ref11;
Pentregue11 = out.potencia11;
ia11 = out.ia11;


figure;
subplot(2,1,1)
plot(tout, w11, 'Color', '#FF0000', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Velocidade angular nas pás [rad/s]')
grid on
hold on
plot(tout, w_ref11, 'g', 'LineWidth', 1.5, 'LineStyle', ':')

subplot(2,1,2)
plot(tout, Ea11, 'Color', '#FF0000', 'LineWidth', 1.5)
xlabel('Tempo [s]')
ylabel('Entradas de controle')
grid on 
hold on
plot(tout, beta11*(180/pi)*100, 'Color', '#FF5C00', 'LineWidth', 1.5)


vecPpe11 = ones(size(tout));
vecPpm11 = ones(size(tout));
vecPg11 = ones(size(tout));



Ppe11 = ((ia11 .* ia11) * Ra); % Perdas elétricas

Ppm11 = (b* (w11 .* w11)); % Perdas mecânicas
 
Pg11 = (kv*ia11.*w11); % Potência gerada

Pe11 = Pg11 - Ppe11 - Ppm11; % Potência entregue
 

figure;
grid on 
hold on
xlabel('Tempo [s]')
ylabel('Potências do sistema [W]')
plot(tout, Ppe11, 'Color', '#FFF700', 'LineWidth', 2.0,'DisplayName', 'Perdas Elétricas')
plot(tout, Ppm11, 'Color', '#FF5C00', 'LineWidth', 2.0, 'DisplayName', 'Perdas Mecânicas')
plot(tout, Pg11, 'Color', '#50C878', 'LineWidth', 2.0, 'DisplayName', 'Potência Gerada')
plot(tout, Pentregue11, 'Color', '#D20A2E', 'LineWidth', 2.0, 'DisplayName', 'Potência Entregue')
legend show;


 


%%









