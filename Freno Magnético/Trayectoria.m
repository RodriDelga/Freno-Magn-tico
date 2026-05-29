function zm = Trayectoria(Bz,z,mag,m,zo,dt,vz_init,gamma,magnet)
    % =========================================================================
    % SIMULACIÓN DE LA TRAYECTORIA DEL IMÁN (CON Y SIN FRENO MAGNÉTICO)
    % =========================================================================
    
    % --- 1. INICIALIZACIÓN DE VARIABLES Y CONDICIONES INICIALES ---
    w = -m * 9.81;      % Cálculo del peso del objeto (Fuerza de gravedad en Newtons)
    zm(1) = zo;         % Vector de posición Z para el imán real (Inicia en zo)
    zmfree(1) = zo;     % Vector de posición Z para la caída libre teórica
    vz(1) = vz_init;    % Vector de velocidad Z para el imán real
    vzfree(1) = 0;      % Velocidad inicial de la caída libre (Parte desde el reposo)
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
                @(z_pos, v_val) a_total_eddy(z_pos, v_val, Bz, z_axis, mag, gamma, m));
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
    
    % --- 3. RECORTE DE DATOS EXCLUSIVO PARA LA CAÍDA LIBRE ---
    % Buscamos los índices numéricos donde el vector de tiempo 'tt' es menor o igual a 1.5s
    idx_recorte = find(tt <= 1.5);
    
    % --- 4. GENERACIÓN DE GRÁFICAS ---
    figure('Color', 'w'); % Creamos una ventana de gráfico con fondo blanco
    
    % Grafica la caída libre RECORTADA (solo los datos correspondientes a los primeros 1.5 segundos)
    plot(tt(idx_recorte), zmfree(idx_recorte), 'r--', 'LineWidth', 2); 
    hold on; % Mantiene la gráfica activa para poder encimar la siguiente línea
    
    % Grafica la trayectoria del imán COMPLETA (hasta que el bucle se detuvo en Z = -3)
    plot(tt, zm, 'g-', 'LineWidth', 2);     
    
    % --- 5. ESTÉTICA Y DETALLES DE LA GRÁFICA ---
    grid on; % Activa la cuadrícula de fondo
    xlabel('Tiempo (s)'); % Etiqueta del eje horizontal
    ylabel('Posición en Z (m)'); % Etiqueta del eje vertical
    title(['Trayectoria del Imán (Magnet = ', num2str(magnet), ')']); % Título dinámico
    legend('Caída libre (0 a 1.5s)', 'Trayectoria magnética completa', 'Location', 'best'); % Leyenda explicativa
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
    k2v = a_func(z + 0.5 * dt * k1v, v + 0.5 * dt * k1v); % Pendiente de velocidad 2: Aceleración a mitad del paso
    
    k3z = v + 0.5 * dt * k2v;          % Pendiente de posición 3: Segunda estimación de velocidad a mitad del paso
    k3v = a_func(z + 0.5 * dt * k2v, v + 0.5 * dt * k2v); % Pendiente de velocidad 3: Segunda aceleración a mitad del paso
    
    k4z = v + dt * k3v;                % Pendiente de posición 4: Velocidad estimada al final del paso dt
    k4v = a_func(z + dt * k3v, v + dt * k3v); % Pendiente de velocidad 4: Aceleración estimada al final del paso dt
    
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

function a = a_total_eddy(z,v,Bz,z_axis,sigma_eff,gamma,m)
    % =========================================================================
    % MODELO DINÁMICO 2: CORRIENTES INDUCIDAS / DE FOUCAULT (Eddy Currents)
    % =========================================================================
    delta = 0.005; % Distancia diferencial para calcular el gradiente espacial
    
    % Interpolación: Misma búsqueda de campo magnético adelante y atrás de la posición
    Bz_foward = interp1(z_axis, Bz, z + delta, "linear", "extrap");
    Bz_backward = interp1(z_axis, Bz, z - delta, "linear", "extrap");
    
    % Cálculo exacto del gradiente del campo magnético (dBz / dz)
    dBz_dz = (Bz_foward - Bz_backward) / (2 * delta);
    
    % CÁLCULO DE LAS FUERZAS FÍSICAS (Ley de Faraday-Lenz)
    % El freno por corrientes inducidas depende del CUADRADO del gradiente del campo y de la velocidad
    F_eddy = -sigma_eff * (dBz_dz^2) * v; 
    Ff = -gamma * v;                    % Fuerza de fricción del aire
    F = F_eddy + Ff - m * 9.81;         % Sumatoria de fuerzas totales: Freno Eddy + Fricción - Peso
    
    a = F / m;                          % Segunda Ley de Newton: Aceleración = Fuerza total / Masa
end