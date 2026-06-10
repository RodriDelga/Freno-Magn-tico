function phiB = flujoB(Bz_slice, x, y, R)
    [X, Y] = meshgrid(x, y);
    mask = (X.^2 + Y.^2) <= R^2;
    mask = mask';                          % alinea con dimensiones de Bz_slice [Lx x Ly]
    dA = abs(x(2)-x(1)) * abs(y(2)-y(1)); % área de cada celda
    phiB = sum(Bz_slice(mask), 'all') * dA;
end


