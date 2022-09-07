%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler
% 1) Resuelve el sistema de tamaño 2 de ODEs
%    y' = f(t,y) en t \in [a,b] sujeta a y(a)=y0
% 2) Calcula errores absolutos.
% 3) Imprime tabla de convergencia en el tiempo en la terminal
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sistema a Resolver
%   f_1(t,y1,y2) = -4*y1  - 2*y2 + cos(t)+ 4*sin(t); 
%   f_2(t,y1,y2) =  3*y1  + y2   - 3*sin(t); 
%   t en el intervalo [0,2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a= 0;
b= 2;
n= 2;        % dimension del sistema de ODEs
N_init= 4;   % numero inicial de subintervalos en el mallado
N_cicles=6;  % numero de ciclos que corre el algoritmo numérico
f_y_t=@f;    %  Declaracion de una funcion
y0= [0, -1]'; %valor inicial @ t=a; vector columna
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Algoritmo
    y_hat=zeros(n,N+1); %Ver Libro de Gilat-Matlab Capitulo 2
                         %El primer indidce es la entrada de la funcion en la
                         %ODE, y el segundo indice es el punto t_i en el intervalo
                         %(a,b)
    y_hat(:,1)=y0; % condicion inicial
        for i=1:N
             y_hat(:,i+1)=y_hat(:,i)+h*f_y_t(t,y_hat(:,i));
            t= t + h;
        end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %** Calcular solución exacta y error utilizando arimética vectorial **
    t_val=linspace(a,b,N+1); %Ver libro Gilat-Matlab pg 38 
    y_val=zeros(n,N+1); %y_val contiene los valores exactos y debe ser de la mismo tamaño que t_hat
    y1_exact_f=@y1_exact; 
    y2_exact_f=@y2_exact; 
    y_val(1,:)=y1_exact_f(t_val);%y_val es un vector renglon de dimension n
    y_val(2,:)=y2_exact_f(t_val);%y_val es un vector renglon de dimension n
    
    %Calcula el error en norma euclideana con respecto a n (tamaño del sistema de ODE)
    e_val= (y_val(1,:)-y_hat(1,:)).^2;         %entrada 1
    e_val= e_val + (y_val(2,:)-y_hat(2,:)).^2; %entrada 2
    e_val=sqrt(e_val);


    %Maximo del error
    err_max(k)=max(e_val); %valor máximo de un vector, ver Gilat-Matlab secc 3.6
    output_text=['El error máximo para N=', num2str(N), ' es: ', num2str(err_max(k))];
    disp(output_text);
end % acaba el ciclo de refinamientos

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Output
% Crear tabla con valores: N| max(y-val - y_hal)| tasa_de_decrecimiento
disp("##### Final del algoritmo #########"); %Imprime en la terminal 
% decrecimiento del error:= log(err_max(i+1)/err_max(i))/log(1/2)
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
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % **La definiciones de las funciones deben ir al final del archivo script
% % Ver Libro de Gilat-Matlab Capitulo 7
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sistema a Resolver
%   f_1(t,y1,y2) = -4*y1  - 2*y2 + cos(t)+ 4*sin(t); 
%   f_2(t,y1,y2) =  3*y1  + y2   - 3*sin(t); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Definicion de la funcion f(t,y)
function [val] = f(t,y)
    val=zeros(2,1); %vector columna
    val(1) = -4*y(1)- 2*y(2)+ cos(t)+ 4*sin(t); 
    val(2) =  3*y(1)+y(2)- 3*sin(t); 
end

%Definicion de la solucion exacta
% 
function [val] = y1_exact(t)
    val = 2*exp(-t)-2*exp(-2*t)+sin(t); 
end


function [val] = y2_exact(t)
    val = -3*exp(-t)+2*exp(-2*t);      
end


