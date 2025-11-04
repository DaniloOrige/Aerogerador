function G = funTeste(w, funTau, b, kv, Ra, Ea)

% função que a raiz retorna w

G = funTau(w) - b*w - (kv/Ra)*(kv*w - Ea);

end


