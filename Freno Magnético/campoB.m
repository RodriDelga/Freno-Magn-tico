function [z] = campoB(ds,km,Px,Py,Pz,dx,dy,nl,N,rw)
    z = -5.2:ds:5.2; % Limites de la malla
    x = z;
    y = z;
    
    Lx = length(x); % Puntos en x
    Ly = length(y); % Puntos en y
    Lz = length(z); % Puntos en z
    % Inicializar los valores en 0 de los diferenciales de las componentes
    % del campo
    dBx = zeros(Lx,Ly,Lz,'single'); 
    dBy = zeros(Lx,Ly,Lz,'single');
    dBz = zeros(Lx,Ly,Lz,'single');
    

    for i = 1:Lx
       for j = 1:Lx
          for k = 1:Lx
             for l = 1:N*nl
                % Componente de r, medidos entre la distancia que hay entre
                % cada punto y el lugar donde se mide
                rx = x(i) - Px(l);
                ry = y(j) - Py(l);
                rz = z(k) - Pz(l);
                % Magnitud de r
                r = sqrt(rx^2 + ry^2 + rz^2 + rw^2);
                r3 = r^3;

                % Ir haciendo la suma de cada componente delñ vector 
                % resultante por medio del producto cruz de dl x r, y 
                % con esas componente, solo multiplicar por la consntante
                dBx(i,j,k) = dBx(i,j,k) + km * dy(l) * rz / r3;
                dBy(i,j,k) = dBy(i,j,k) + km * dx(l) * rz / r3;
                dBz(i,j,k) = dBz(i,j,k) + km * (dx(l) * ry - dy(l) * rx) / r3;
             end
          end
       end 
    end

% Corte en y = 0 → plano xz

% Vector de coordenadas del grid (el mismo que usa campoB internamente)
coords = -5.2:ds:5.2;

% Encuentra el índice central del vector, que corresponde a y ≈ 0
% round() por si length es par y no cae exacto en el centro
iy = round(length(coords) / 2);

% Crea una malla 2D en el plano xz
% X2 varía en columnas (dirección x), Z2 varía en filas (dirección z)
% Ambas matrices son de tamaño length(coords) x length(coords)
[X2, Z2] = meshgrid(coords, coords);

% Extrae el plano xz del array 3D dBx fijando el índice de y en iy
% squeeze() elimina la dimensión colapsada (que queda de tamaño 1)
% El resultado sin transponer sería (Lx × Lz), pero meshgrid espera (Lz × Lx)
% por eso se transpone con '
Bx_slice = squeeze(dBx(:, iy, :))';
Bz_slice = squeeze(dBz(:, iy, :))';

% Calcula la magnitud del campo en cada punto del plano
Bmag = sqrt(Bx_slice.^2 + Bz_slice.^2);

% Figura

% Crea una figura con fondo blanco
figure('Color', 'w')

% Dibuja un mapa de color relleno con 40 niveles de contorno
% 'LineColor','none' quita las líneas entre contornos para que se vea suave
contourf(X2, Z2, Bmag, 40, 'LineColor', 'none')

% Paleta de colores: negro→azul→verde→amarillo→rojo (buena para campos)
colormap turbo

% Agrega la barra de color y le pone etiqueta
cb = colorbar;
cb.Label.String = '|B|';

% Permite agregar más elementos encima sin borrar el contourf
hold on

% Dibuja las líneas de campo a partir de los vectores (Bx, Bz)
% El último argumento (2) controla la densidad de líneas: mayor = más líneas
streamslice(X2, Z2, Bx_slice, Bz_slice, 2)

% Etiquetas de ejes y título
xlabel('x')
ylabel('z')
title('Campo B — plano xz (y = 0)')

% axis equal → misma escala en x y z para no distorsionar el campo
% tight → recorta los márgenes vacíos alrededor de la gráfica
axis equal tight