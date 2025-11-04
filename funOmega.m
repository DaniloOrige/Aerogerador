function dwdt = funOmega(J, tau, b, w, kv, ia)

dwdt = (1/J)*(tau - b*w - kv*ia);

end