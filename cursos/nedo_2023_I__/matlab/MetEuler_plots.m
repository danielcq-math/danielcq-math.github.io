%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler
% 1) Resuelve la ODE
%    y' = f(t,y) en t \in [a,b] sujeta a y(a)=y0
% 2) Calcula errores relativos y absolutos.
% 3) Grafica la solucion exacta y aproximada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=2;
N=4; %#total de puntos en el mallado
h=(b-a)/N;
y0= 0.5; %valor inicial @ t=a
f_y_t=@f; %Declaracion de una funcion
t=a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Algoritmo
y_hat=zeros(N+1,1); %Ver Libro de Gilat-Matlab Capitulo 2
y_hat(1)=y0; %El primer elemento de un vector en Matlab de empieza con índice 1 
for i=1:N
    y_hat(i+1)=y_hat(i)+h*f_y_t(t,y_hat(i));
    t= t + h;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcular solución exacta y error
y_val=zeros(N+1,1);
e_val =zeros(N+1,1);
e_val_r=zeros(N+1,1);
t_val=zeros(N+1,1);
t_val(1)=a;
y_exact_f=@y_exact;
y_val(1)=y_exact_f(t_val(1)); %El primer elemento de un vector en Matlab de empieza con índice 1 
for i=1:N
    t_val(i+1)= t_val(i) + h;
    y_val(i+1)=y_exact_f(t_val(i+1));
    e_val(i+1)= abs(y_val(i+1)-y_hat(i+1));
    e_val_r(i+1)=e_val(i+1)/abs(y_val(i+1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(Para más info de como graficar en Matlab ver libro de Gilat Cap 5) 
%Grafica y_hat vs y_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bloque hold on - hold off (Ver Libro Gilat-Matlab Seccion 5.3.2)
hold on %Instrucciones para un grafica
plot(t_val,y_hat,'-b*',t_val,y_val,'-r');% 2 funciones
title('SolApprox vs SolExacta');
legend('y\_hat','y\_val'); % etiquetas de lo que se esta graficando
                           % cuando se utilizan carácteres especiales hay
                           % utilizar primero el caracter '\'

hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure %Crear otra ventana para graficar
 %Grafica del error absoluto
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Bloque hold on - hold off
 hold on
 title('Error Abs');
 plot(t_val,e_val,'-r*')%error absoluto
 legend('e\_val');
 hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure %Crear otra ventana para graficar
 %Grafica del error relativo
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Bloque hold on - hold off
 hold on
 title('Error Rel');
 plot(t_val,e_val_r,'-r*')%error relativo
 legend('e_r\_val');
 hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capítulo 7


%Definición de la función f(t,y)
function [val] = f(t,y)
    val = y-t^2+1; 
end

%Definición de la solución exacta
function [val] = y_exact(t)
    [val] = (t+1)^2 - 0.5*exp(t); 
end
