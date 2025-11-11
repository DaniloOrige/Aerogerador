function F = funNominal(x, rho, A, b, kv, R, Ra)

U = x(1);
w = x(2);
Ea = x(3);
ia = x(4);

Cpmax = 0.44;
Pw = 0.5*rho*A*U^3;
Pnom = 4.2*10^6;
lambdaopt = 7.5;

eq1 = (Cpmax * Pw)/w - b*w - kv*ia; 
eq2 = lambdaopt - (w*R)/U;
eq3 = kv*w - Ea - Ra*ia;
eq4 = Pnom - Ea*ia;


F = [eq1; eq2; eq3; eq4];



