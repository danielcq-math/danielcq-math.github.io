% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler Implicito
% 1) Resuelve la ODE
%    y' = f(t,y) en t \in [a,b] sujeta a y(a)=y0
% 2) Calcula errores relativos y absolutos.
% **Utiliza aritmética vectorial**
% Imprime tabla de convergencia en el tiempo en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a=0;
b=1;
N_init=3;          % numero inicial del total de sub-intervalos en el mallado
N_cicles=8;         % numero de ciclos que corre el algoritmo numérico
f_y_t=@f; %Declaracion de una funcion
y0= -1; %valor inicial @ t=a
%Parametros PFijo
N_iter_PFijo=100; %Numero de iteraciones del Algoritmo del PF
eps_tol= 1.e-10;  %Tol abosluto del criterio de paro del Alg. del PF
format shortEng

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construimos un vector que obtiene el numero total de puntos en el mallado
%para cada ciclo de refinamiento 

N_vec=zeros(1,N_cicles);
N_vec(1)=N_init; %El primer ciclo no se refina

for k=2:N_cicles
  N_vec(k)=2*N_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end

y_hat_max=zeros(1,N_cicles); % Vector que contiene el maximo del valor absoulto de y_hat para cada ciclo
err_max=zeros(1,N_cicles); % Vector que contiene el maximo del error para cada ciclo
% ciclos de refinamiento
for k=1:N_cicles
    N=N_vec(k); %actualizando el numero de subintervalos
    h=(b-a)/N;
    t=a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Algoritmo
    y_hat=zeros(N+1,1); %Ver Libro de Gilat-Matlab Capitulo 2   
    y_hat(1)=y0; %El primer elemento de un vector en Matlab de empieza con índice 1 
    for i=1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Metodo del Punto fijo para Euler implicito
        p0=y_hat(i); %Aproximacion inicial
        j=1;
        while j<N_iter_PFijo    
            p=y_hat(i) + h*f_y_t(t+h,p0);% t+h corresponde al valor de y_hat(i+1)
            err_p_fijo = abs(p-p0);
            p0=p;
            %disp(err_p_fijo);
            if(err_p_fijo < eps) %Criterio de paro
                break %sale del ciclo while
            end  
            j=j+1;
        end
        if( j==N_iter_PFijo)
            disp("** Algoritmo de PFijo falló**")
            break
        end
        y_hat(i+1)=p;
        t= t + h;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%** Calcular solución exacta y error utilizando arimética vectorial **
    t_val=linspace(a,b,N+1); %Ver libro Gilat-Matlab pg 38 
    y_exact_f=@y_exact;
    y_val=y_exact_f(t_val);
    e_val= abs(y_val-transpose(y_hat)); %y_val es un vector-fila y y_hat un vector columna

    y_hat_max(k)=max(abs(y_hat)); %valor máximo
    output_text=['El valor máximo de y_hat para N=', num2str(N), ' es: ', num2str(y_hat_max(k))];
    disp(output_text);
    err_max(k)=max(e_val); %valor máximo de un vector, ver Gilat-Matlab secc 3.6
    output_text=['El error máximo para N=', num2str(N), ' es: ', num2str(err_max(k))];
    disp(output_text);


end %acaba el ciclo de refinamientos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% Crear tabla con valores: N| max(y-val - y_hal)| tasa_de_decrecimiento
disp("##### Final del algoritmo #########"); %Imprime en la terminal 
% decrecimiento del error:= log(err_max(i+1)/err_max(i))*log(2)
err_rate=zeros(1,N_cicles); %tasa de decrecimiento del error
err_rate(1)=1;
for i=2:N_cicles
    err_rate(i)=log(err_max(i)/err_max(i-1))/log(1/2);
end
output_table=[N_vec' y_hat_max' err_max' err_rate'] ; % Crea tabla con vectores columna
% %Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla: N_vec'  y_hat_max' err_max' err_rate'"); %Imprime en la terminal                                                 
disp(output_table);%Imprime en la terminal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7


%Definicion de la funcion f(t,y)
function [val] = f(t,y)
    val = 7*exp(7*t).*(y-0.25*t).^2 + 0.25; 
end

%Definicion de la solucion exacta
function [val] = y_exact(t)
    [val] = 0.25*t - exp(-7*t); 
end
