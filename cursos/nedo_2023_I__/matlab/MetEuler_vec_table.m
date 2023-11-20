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
% **Utiliza aritmética vectorial**
% Imprime un tabla en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=2;
N=10; %#total de puntos en el mallado
f_y_t=@f; %Declaracion de una funcion
y0= 0.5; %valor inicial @ t=a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(b-a)/N;
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
%** Calcular solución exacta y error utilizando arimética vectorial **
t_val=linspace(a,b,N+1); %Ver libro Gilat-Matlab pg 38 
y_exact_f=@y_exact;
y_val=y_exact_f(t_val);
e_val= abs(y_val-transpose(y_hat)); %y_val es un vector-fila y y_hat un vector columna
e_val_r=e_val./abs((y_val));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% Crear tabla con valores: t_val| y_hal| y_val| errores abs | errores rel
% Crea tabla con vectores columna
output_table=[t_val' y_hat y_val' e_val' e_val_r'] ; % y_hat es el unico vector columna no necesita tranpuesta
%Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla: t_val | y_hal|  y_exact| errores abs | errores rel:"); %Imprime en la terminal                                                 
disp(output_table);%Imprime en la terminal


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
