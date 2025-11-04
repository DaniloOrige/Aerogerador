function diadt = funIa(La, kv, w, Ra, Ea, ia)

diadt = (1/La)*(kv*w - Ra*ia - Ea);

end