%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler Explicito (exploracion de estabilidad)
% Resuelve la ODE
%    y' = lambda*y en t \in [a,b] sujeta a y(a)=y0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=1;
N=10; %#total de intervalos que subdivide a [a,b] 
y0= 1; %valor inicial @ t=a
lambda=-50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Algoritmo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inicializacion
t=a;
delta_t=(b-a)/N;
y_hat=zeros(N+1,1); %Ver Libro de Gilat-Matlab Capitulo 2
y_hat(1)=y0; %El primer elemento de un vector en Matlab de empieza con índice 1 
%Calcular solución aproximada
for i=1:N  %#total de intervalos que subdivide a [a,b]
    y_hat(i+1)=y_hat(i)+delta_t*lambda*y_hat(i);
    t= t + delta_t;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica del solucion aproximada y la analitica
y_exact_f=@y_exact; %sol analitica
t_vec=linspace(a,b,N+1);
y_exact_vec=y_exact_f(t_vec,lambda);
plot(t_vec,y_hat,'*-b',t_vec,y_exact_vec,'-r');
%plot(t_vec,y_exact_vec,'-r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7


%Definicion de la solucion exacta
function [val] = y_exact(t,lambda)
    [val] = exp(lambda*t); 
end