function [dBz_eje, z,dPhi_dz] = campoB(ds,km,Px,Py,Pz,dx,dy,nl,N,rw,plot_option,R)
    % 1. Definición de la malla según la opción elegida
    z = -5.2:ds:5.2; % Límites de la malla en Z (siempre fijos para la trayectoria)
    
    if plot_option
        x = z;
        y = z;
    else
        % Malla mínima en X e Y para ahorrar memoria y tiempo de CPU
        x = -R:0.1:R;
        y = -R:0.1:R;
    end
    
    Lx = length(x); 
    Ly = length(y); 
    Lz = length(z); 
    
    % Inicializar matrices en 0
    dBx = zeros(Lx,Ly,Lz,'single'); 
    dBy = zeros(Lx,Ly,Lz,'single');
    dBz_malla = zeros(Lx,Ly,Lz,'single');
    
    % --- CORRECCIÓN ABSOLUTA DE LOS BUCLES ---
    for i = 1:Lx
       for j = 1:Ly  % Corregido a Ly
          for k = 1:Lz  % Corregido a Lz
             for l = 1:N*nl
                rx = x(i) - Px(l);
                ry = y(j) - Py(l);
                rz = z(k) - Pz(l);
                
                r = sqrt(rx^2 + ry^2 + rz^2 + rw^2);
                r3 = r^3;
                
                dBx(i,j,k) = dBx(i,j,k) + km * dy(l) * rz / r3;
                dBy(i,j,k) = dBy(i,j,k) + km * dx(l) * rz / r3;
                dBz_malla(i,j,k) = dBz_malla(i,j,k) + km * (dx(l) * ry - dy(l) * rx) / r3;
             end
          end
       end 
    end


    % 2. Módulo de graficación del plano XZ (Solo si plot_option es verdadero)
    if plot_option
        coords = -5.2:ds:5.2;
        iy = round(length(coords) / 2);
        [X2, Z2] = meshgrid(coords, coords);
        
        Bx_slice = squeeze(dBx(:, iy, :))';
        Bz_slice = squeeze(dBz_malla(:, iy, :))';
        Bmag = sqrt(Bx_slice.^2 + Bz_slice.^2);
        
        figure('Color', 'w', 'Name', 'Plano XZ')
        contourf(X2, Z2, Bmag, 40, 'LineColor', 'none')
        colormap turbo; cb = colorbar; cb.Label.String = '|B|';
        hold on
        streamslice(X2, Z2, Bx_slice, Bz_slice, 2)
        xlabel('x (m)'); ylabel('z (m)');
        title('Campo B — plano xz (y = 0)');
        axis equal tight
    end

    % 3. Extracción del perfil central (Eje de simetría X=0, Y=0)
    idx_x = ceil(Lx / 2);
    idx_y = ceil(Ly / 2);
    dBz_eje = squeeze(dBz_malla(idx_x, idx_y, :)); % Vector 1D que va al Main
    
    for k = 1:length(z)
        Bz_slice = squeeze(dBz_malla(:,:,k));
        phiB(k) = flujoB(Bz_slice,x,y,R);
    end

    dPhi_dz = diff(phiB) ./ diff(z);

    % 4. Gráfica del Gradiente del Campo (dBz/dz) contra Z
    % Usamos gradient para no perder dimensiones ni desfasar el vector z
    dz_espacio = z(2) - z(1);
    dBz_dz_profile = gradient(dBz_eje, dz_espacio);
    
    figure('Color', 'w', 'Name', 'Gradiente Bz')
    plot(z, dBz_dz_profile, 'r-', 'LineWidth', 2);
    xlabel('z (m)'); ylabel('dB_z / dz (T/m)');
    title('Gradiente del campo magnético B_z');
    grid on;
end