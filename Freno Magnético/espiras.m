function [Px,Py,Pz,dx,dy,dz] = espiras(sz,nl,N,R)
    dtheta = 2*pi / N; % Distancia entre los puntos
    ang = 0:dtheta:(2*pi-dtheta); % Vector con los angulos

    s = 1; % variable para cambiar los puntos de s
    for i = 1:nl
        Px(s:s+N-1) = R * cos(ang); % Componentes en x
        Py(s:s+N-1) = R * sin(ang); % Componentes en y
        Pz(s:s+N-1) = -(nl - 1)/2*sz + (i-1)*sz; % Componentes en z posicionadas

        dx(s:s+N-1) = -Py(s:s+N-1) * dtheta; % Calcular los diferenciales en x
        dy(s:s+N-1) =  Px(s:s+N-1) * dtheta; % Calcular los diferenciales en y
        s = s + N; % Actualizar el valor de s
    end
    dz = zeros(1, N*nl); % Initialize z-direction displacement
    figure(1)
    quiver3(Px,Py,Pz,dx,dy,dz,0.5,'-r','LineWidth',2)
    view(-34,33)
end