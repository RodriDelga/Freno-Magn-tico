addpath('C:\Users\yared\Documents\MATLAB\Freno Magnético')

nl = 5; % número de espiras
N = 20; % Puntos en las espiras
R = 1.5; % Radio
sz = 1; % Distancia entre espiras
I = 300; % Corriente
mo = 4*pi*1e-7; % Permeanilidad magnética
km = mo * I / (4*pi); % Constante de Biot-Savart
rw = 0.2; % Constante de Biot-Savart
ds = 0.1; % Grosor entre los puntos de la malla

[Px,Py,Pz,dx,dy,dz] = espiras(sz,nl,N,R);
[z] = campoB(ds,km,Px,Py,Pz,dx,dy,nl,N,rw);