function zm = Trayectoria(Bz,z,mag,m,zo,dt,vz_init,gamma,magnet,dPhi_dz,R_circuito)
    % =========================================================================
    % SIMULACIÓN DE LA TRAYECTORIA DEL IMÁN (CON Y SIN FRENO MAGNÉTICO)
    % =========================================================================
    
    % --- 1. INICIALIZACIÓN DE VARIABLES Y CONDICIONES INICIALES ---
    w = -m * 9.81;      % Cálculo del peso del objeto (Fuerza de gravedad en Newtons)
    zm(1) = zo;         % Vector de posición Z para el imán real (Inicia en zo)
    zmfree(1) = zo;     % Vector de posición Z para la caída libre teórica
    vz(1) = vz_init;    % Vector de velocidad Z para el imán real
    vzfree(1) = vz_init;      % Velocidad inicial de la caída libre (Parte desde el reposo)
    tt(1) = 0;          % Vector de tiempo (Arranca en 0 segundos)
    cc = 1;             % Contador/Índice que controla los pasos dentro del bucle
    z_axis = z;         % Guardamos el vector de posiciones de la malla del campo magnético
    
    % --- 2. BUCLE PRINCIPAL DE SIMULACIÓN (CONTROLADO POR ALTURA) ---
    % El bucle se ejecutará continuamente mientras el imán real esté arriba de Z = -3 metros
    while zm(cc) > -3
        
        % A) CAÍDA LIBRE TEÓRICA (Método de Euler estándar)
        % Nueva posición en base a la velocidad actual y la aceleración de la gravedad
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*(-9.81)*dt^2;
        % Nueva velocidad bajo la aceleración constante de la gravedad (-9.81 m/s^2)
        vzfree(cc+1) = vzfree(cc) - 9.81 * dt;
        
        % B) CAÍDA DEL IMÁN REAL (Cálculo mediante Runge-Kutta 4)
        if magnet
            % Si magnet es verdadero (1), calcula usando la interacción dipolar clásica
            [zm(cc+1), vz(cc+1)] = rk4_step(zm(cc), vz(cc), dt, ...
                @(z_pos, v_val) a_total(z_pos, v_val, Bz, z_axis, mag, gamma, m));
        else
            % Si magnet es falso (0), calcula usando el modelo de corrientes de Foucault (Eddy)
            [zm(cc+1), vz(cc+1)] = rk4_step(zm(cc), vz(cc), dt, ...
                @(z_pos, v_val) a_total_eddy(z_pos,v_val,Bz,z_axis,gamma,m,dPhi_dz,R_circuito));
        end
        
        % C) ACTUALIZACIÓN DE TIEMPO Y CONTADORES
        tt(cc+1) = tt(cc) + dt;  % Avanzamos el reloj un diferencial de tiempo 'dt'
        cc = cc + 1;             % Incrementamos el índice para calcular el siguiente paso
        
        % D) CONDICIÓN DE PARADA PREMATURA
        % Si el imán cae tan lento que su velocidad es casi cero (menor a 1 mm/s), detiene el bucle
        if abs(vz(cc)) < 1e-3
            disp("velocidad pequeña")
            break; % Rompe el bucle 'while' de forma inmediata
        end
    end
    
    vz_out = vz;

    % --- RECORTE CAÍDA LIBRE ---
    idx_recorte = find(tt <= 1.5);

    % --- GRÁFICA POSICIÓN ---
    figure('Color', 'w');
    plot(tt(idx_recorte), zmfree(idx_recorte), 'r--', 'LineWidth', 2);
    hold on;
    plot(tt, zm, 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Tiempo (s)');
    ylabel('Posición en Z (m)');
    title(['Trayectoria del Imán (Magnet = ', num2str(magnet), ')']);
    legend('Caída libre (0 a 1.5s)', 'Trayectoria magnética completa', 'Location', 'best');

    % --- GRÁFICA VELOCIDAD ---
    figure('Color', 'w');
    plot(tt, vz, 'g-', 'LineWidth', 2);
    grid on;
    hold on;
    plot(tt(idx_recorte), vzfree(idx_recorte), 'b--', 'LineWidth', 2);
    xlabel('Tiempo (s)');
    ylabel('Velocidad en Z (m/s)');
    title(['Velocidad del Imán (Magnet = ', num2str(magnet), ')']);
    legend('Velocidad magnética', 'Location', 'best');
end

function [z_next, v_next] = rk4_step(z,v,dt,a_func)
    % =========================================================================
    % SOLVER NUMÉRICO: RUNGE-KUTTA DE 4º ORDEN (RK4)
    % =========================================================================
    % Calcula 4 pendientes distintas (k1 a k4) en intervalos de tiempo intermedios
    % para ofrecer una aproximación de la posición y velocidad extremadamente precisa.
    
    k1z = v;                           % Pendiente de posición 1: Velocidad actual
    k1v = a_func(z,v);                 % Pendiente de velocidad 1: Aceleración actual
    
    k2z = v + 0.5 * dt * k1v;          % Pendiente de posición 2: Velocidad estimada a mitad del paso
    k2v = a_func(z + 0.5 * dt * k1z, v + 0.5 * dt * k1v); % Pendiente de velocidad 2: Aceleración a mitad del paso
    
    k3z = v + 0.5 * dt * k2v;          % Pendiente de posición 3: Segunda estimación de velocidad a mitad del paso
    k3v = a_func(z + 0.5 * dt * k2z, v + 0.5 * dt * k2v); % Pendiente de velocidad 3: Segunda aceleración a mitad del paso
    
    k4z = v + dt * k3v;                % Pendiente de posición 4: Velocidad estimada al final del paso dt
    k4v = a_func(z + dt * k3z, v + dt * k3v); % Pendiente de velocidad 4: Aceleración estimada al final del paso dt
    
    % FÓRMULA DE ACTUALIZACIÓN DE RK4: Combinación ponderada de las 4 pendientes
    z_next = z + dt/6 * (k1z + 2*k2z + 2*k3z + k4z); % Siguiente posición calculada
    v_next = v + dt/6 * (k1v + 2*k2v + 2*k3v + k4v); % Siguiente velocidad calculada
end

function a = a_total(z,v,Bz,z_axis,mag,gamma,m)
    % =========================================================================
    % MODELO DINÁMICO 1: ACOPLAMIENTO DIPOLAR CLÁSICO (Fuerza del imán)
    % =========================================================================
    delta = 0.005; % Distancia diferencial para calcular el gradiente espacial
    
    % Interpolación: Buscamos cuánto vale el campo magnético un paso adelante y un paso atrás del imán
    Bz_foward = interp1(z_axis, Bz, z + delta, "linear", "extrap");
    Bz_backward = interp1(z_axis, Bz, z - delta, "linear", "extrap");
    
    % Diferencias finitas centrales para obtener el gradiente del campo magnético (dBz / dz)
    dBz_dz = (Bz_foward - Bz_backward) / (2 * delta);
    
    % CÁLCULO DE LAS FUERZAS FÍSICAS
    Fm = -mag * dBz_dz;       % Fuerza magnética dipolar (F = -m * dB/dz)
    Ff = -gamma * v;          % Fuerza de fricción viscosa del aire (Opuesta a la dirección del movimiento)
    F = Fm + Ff - m * 9.81;   % Sumatoria de fuerzas totales: Magnética + Fricción - Peso
    
    a = F / m;                % Segunda Ley de Newton: Aceleración = Fuerza total / Masa
end

function a = a_total_eddy(z,v,Bz,z_axis,gamma,m,dPhi_dz,R_circuito)

    z_mid = 0.5 * (z_axis(1:end-1) + z_axis(2:end));
    
    % Forzamos que z sea escalar y el resultado también lo sea
    dPhi_dz_local = interp1(z_mid(:), dPhi_dz(:), z, 'linear', 0);
    dPhi_dz_local = dPhi_dz_local(1);  % garantiza escalar estricto
    
    fem     = -dPhi_dz_local * v;
    I_ind   = fem / R_circuito;
    F_eddy  = I_ind .* dPhi_dz_local;  % .* por si acaso
    Ff      = -gamma * v;
    Fg      = -m * 9.81;
    
    a = (F_eddy + Ff + Fg) / m;
end

