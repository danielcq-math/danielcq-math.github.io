%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler
% Resuelve la ODE
%    y' = f(t,y) en t \in [a,b] sujeta a y(a)=y0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=2;
N=10; %#total de puntos en el mallado
y0= 0.5; %valor inicial @ t=a
f_y_t=@f; %Declaracion de una funcion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Algoritmo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inicializacion
t=a;
h=(b-a)/N;
y_hat=zeros(N+1,1); %Ver Libro de Gilat-Matlab Capitulo 2
y_hat(1)=y0; %El primer elemento de un vector en Matlab de empieza con índice 1 
%Calcular solución aproximada
for i=1:N
    y_hat(i+1)=y_hat(i)+h*f_y_t(t,y_hat(i));
    t= t + h;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7

%Definicion de la funcion f(t,y)
function [val] = f(t,y)
    val = y-t^2+1; 
end

%Definicion de la solucion exacta
function [val] = y_exact(t)
    [val] = (t+1)^2 - 0.5*exp(t); 
end
