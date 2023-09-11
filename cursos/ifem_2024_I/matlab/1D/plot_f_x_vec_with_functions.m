%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso Introduccion a los Elementos Finitos
% Lic. en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
% utilizando aritmetica vectorial
% y definiendo funciones en matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input, parametros del usuario
N= 4; %Numero de subintervalos que dividen a [a,b]
a=0;
b=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
f_x=@function_y; %Libro de Gilat-Matlab Capitulo 7
g_x=@function_z; 
x=linspace(a,b,N+1);%Crea un vector FILA de N+1 puntos que van de [a,b]
                  % y cada punto esta separado por h=(b-a)/N
                  % Ver libro Gilat-Matlab pg 38 
y=f_x(x);  
z=g_x(x);
plot(x,y,'-r*',x,z, '-b*'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
%   Ver Libro de Gilat-Matlab Capitulo 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funciones definidas por el usuario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: function_y
%       Input:
%              x: vector de valores (Dominio)  

% Descripcion:
%             Funcion a graficar
function [val] = function_y(x)
                 val=x.^2; %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: function_z
%       Input:
%              x: vector de valores (Dominio)  

% Descripcion:
%             Funcion a graficar
function [val] = function_z(x)
                 val=x.^3; %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
end