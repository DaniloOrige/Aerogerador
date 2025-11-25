function [dxdt] = funDin_aerogerador(J, U_nom, Cpmax, rho, A, b, La, Ra, Ea, kv, x)

ia = x(1);
w = x(2);

Pw = 0.5*rho*A*U_nom^3;

diadt = (1/La) * (kv*w - Ra*ia - Ea);

dwdt = (1/J) * ((Cpmax * Pw)/w - b*w - kv*ia);


dxdt = [diadt; dwdt];

end

