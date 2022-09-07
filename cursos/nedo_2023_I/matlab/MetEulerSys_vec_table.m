%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Método de Euler
% 1) Resuelve el sistema de tamaño 2 de ODEs
%    y_' = f_(t,y) en t \in [a,b] sujeta a y_(a)=y_0
% 2) Calcula errores relativos y absolutos.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sistema a Resolver
%   f_1(t,y1,y2) = -4*y1  - 2*y2 + cos(t)+ 4*sin(t); 
%   f_2(t,y1,y2) =  3*y1  + y2   - 3*sin(t); 
%   t en el intervalo [0,2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **Utiliza aritmética vectorial**
% Imprime un tabla en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parametros del metodo de Euler
%interval endpoints of [a,b]
a= 0;
b= 2;
n= 2;  %dimension del sistema de ODEs
N= 10; %#total de subintervalos en el mallado
f_y_t=@f; %Declaracion de una funcion
y0= [0, -1]'; %valor inicial @ t=a; vector columna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=(b-a)/N;
t=a;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Algoritmo
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

e_val= (y_val(1,:)-y_hat(1,:)).^2; %y_val es un vector-fila y y_hat un vector columna
e_val= e_val + (y_val(2,:)-y_hat(2,:)).^2; %y_val es un vector-fila y y_hat un vector columna
e_val=sqrt(e_val);
% sqrt(e_1(ti)*e_1(ti) + e_2(ti)*e_2(ti)) = l2(e_(ti))
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Output
% % Crear tabla con valores: t_val| errores abs 
% % Crea tabla con vectores columna
output_table=[t_val' e_val'] ; % y_hat es el unico vector columna no necesita tranpuesta
% %Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla: t_val | errores abs:"); %Imprime en la terminal                                                 
disp(output_table);%Imprime en la terminal
% 
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


