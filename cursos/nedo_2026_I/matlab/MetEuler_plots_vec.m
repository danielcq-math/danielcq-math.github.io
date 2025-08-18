%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Analisis Numérica de Ecuaciones Diferenciales Ordinarias
% IIMAS-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler
% 1) Resuelve la ODE
%    y' = f(t,y) en t \in [a,b] sujeta a y(a)=y0
% 2) Calcula errores relativos y absolutos.
% 3) Grafica la solucion exacta y aproximada
% **Utiliza aritmética vectorial**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros de usuario del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=2;
N=4; %#total de sub-intervalos en el mallado
y0= 0.5; %valor inicial @ t=a
f_y_t=@f; %Declaracion de una funcion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
t=a;
dt=(b-a)/N;

%Algoritmo
y_hat=zeros(N+1,1); %Ver Libro de Gilat-Matlab Capitulo 2
y_hat(1)=y0; %El primer elemento de un vector en Matlab de empieza con índice 1 
for i=1:N
    y_hat(i+1)=y_hat(i)+dt*f_y_t(t,y_hat(i));
    t= t + dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-procesamiento
%** Calcular solución exacta y error utilizando arimética vectorial **

t_val=linspace(a,b,N+1);%Coinciden con los puntos de la malla, %Ver libro Gilat-Matlab pg 38 
y_exact_f=@y_exact;
y_val=y_exact_f(t_val); % valores de la sol. exacta
e_val= abs(y_val-transpose(y_hat)); %y_val es un vector-fila y y_hat un vector columna
e_val_r=e_val./abs((y_val));



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %(Para más info de como graficar en Matlab ver libro de Gilat Cap 5) 
% %Grafica y_hat vs y_val
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Bloque hold on - hold off (Ver Libro Gilat-Matlab Seccion 5.3.2)
hold on %Instrucciones para un grafica
plot(t_val,y_hat,'-b*',t_val,y_val,'-r');
title('SolApprox vs SolExacta');
legend('y\_hat','y\_val');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %Crear otra ventana para graficar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica del error absoluto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bloque hold on - hold off
hold on
title('Error Abs');
plot(t_val,e_val,'-r*')
legend('e\_val');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7


%Definicion de la funcion f(t,y)
function [val] = f(t,y)
    val = y-t.^2+1; 
end

%Definicion de la solucion exacta
function [val] = y_exact(t)
    [val] = (t+1).^2 - 0.5*exp(t); 
end
