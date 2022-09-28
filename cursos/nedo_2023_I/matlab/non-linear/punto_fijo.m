%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.
% castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ilustracion del método del punto fijo
%
% Buscamos la solución única de la ecuación
%  x^3 + 4x^2 -10 = 0 en el intervalo [1,2]
% Para encontrar esta solucion utilizamos el algortimo del punto fijo
% para 
%       g1(x)= x+x^3+4x^2-10
%       g2(x)= sqrt(10/x -4x)
%       g3(x)= 0.5*sqrt(10-x^3) pero con p0=g3(1.5)
%       g4(x)= sqrt(10/4+x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros
N0=100; %numero máximo de iteraciones
g_x=@g4; %Declaracion de una funcion
a=1;
b=2;
%p0=g_x(a); %p0 es la aproximacion inicial del punto fijo
%p0=g_x(b);
p0=g_x((a+b)*0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algoritmo
%Para saber más del bucle while: Ver Gilat-Matlab Sec 6.4.2
i=1;
format longEng %imprimir varios 15-digitos de precision
%format shortEng %imprimir 4-digitos de precision
while i<=N0
    p=g_x(p0); %Calculo de la Sol aproximada
    disp(p); %Imprimir
    p0=p; %Update de la aproximacion
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7
%Definicion de la funcion g1(x)
function [val] = g1(x)
    val = x+x^3+4*x^2-10; 
end

function [val] = g2(x)
    val = x+x^3+4*x^2-10; 
end

function [val] = g3(x)
    val = 0.5*sqrt(10-x^3); 
end

function [val] = g4(x)
    val = sqrt(10/(4+x)); 
end

