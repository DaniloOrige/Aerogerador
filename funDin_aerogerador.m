function [dxdt] = funDin_aerogerador(J, funTau, b, La, Ra, Ea, kv, x)

ia = x(1);
w = x(2);

tau = funTau(w);

diadt = (1/La) * (kv*w - Ra*ia - Ea);

dwdt = (1/J) * (tau - b*w - kv*ia);


dxdt = [diadt; dwdt];

end

