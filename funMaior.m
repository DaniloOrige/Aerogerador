function beta = funMaior (lambdaopt, U, R, rho, A, ia, kv, beta, b)

Pw = 0.5*rho*A*U^3;

w = (U*lambdaopt)/R;
Cp = (0.44 - 0.167*beta)*sin((pi*lambdaopt)/(15 - 0.3*beta)) - 0.16*lambdaopt*beta;
beta = (Cp * Pw)/w - b*w - kv*ia; 

end