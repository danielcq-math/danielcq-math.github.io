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
% Imprime tabla de convergencia en el tiempo en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=2;
c=1;             % Cálculo de la solucion  approx en t=c
N_init=2;        % numero inicial del total de sub-intervalos en el mallado
N_cicles=8;      % numero de ciclos que corre el algoritmo numerico
f_y_t=@f; %Declaracion de una funcion
y0= 0.5; %valor inicial @ t=a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construimos un vector que obtiene el numero total de puntos en el mallado
%para cada ciclo de refinamiento 
N_vec=zeros(1,N_cicles);
N_vec(1)=N_init; %El primer ciclo no se refina

for k=2:N_cicles
  N_vec(k)=2*N_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end

err_max=zeros(1,N_cicles); % Vector que contiene el maximo del error para cada ciclo
% ciclos de refinamiento
for k=1:N_cicles
    output_text=['## Ciclo de refinamiento k=', num2str(k-1)];
    disp(output_text);

    N=N_vec(k);
    h=(b-a)/N;
    t=a;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    err_max(k)=max(e_val); %valor máximo de un vector, ver Gilat-Matlab secc 3.6
    output_text=['El error máximo para N=', num2str(N), ' es: ', num2str(err_max(k))];
    disp(output_text);
    
    %Aproximar el valor de y(c) usando y_hat por interpolacion
    y_hat_c=interp1(t_val,y_hat,c,'linear'); %ver Gilat-Matlab secc 8.3
    output_text=['El valor approximado de y(c) en c=', num2str(c), ' es: ', num2str(y_hat_c)];
    disp(output_text);
end % acaba el ciclo de refinamientos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Output
% Crear tabla con valores: N| max(y-val - y_hal)| tasa_de_decrecimiento
disp("##### Final del algoritmo #########"); %Imprime en la terminal 
% decrecimiento del error:= log(err_max(i+1)/err_max(i))*log(2)
err_rate=zeros(1,N_cicles); %tasa de decrecimiento del error
err_rate(1)=1;
for i=2:N_cicles
    err_rate(i)=log(err_max(i)/err_max(i-1))/log(1/2);
end
output_table=[N_vec' err_max' err_rate'] ; % Crea tabla con vectores columna
% %Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla: N_vec' err_max' err_rate'"); %Imprime en la terminal                                                 
disp(output_table);%Imprime en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fin del programa
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
