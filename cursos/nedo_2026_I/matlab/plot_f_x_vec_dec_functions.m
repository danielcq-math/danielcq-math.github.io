%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
% utilizando aritmetica vectorial de matlab y funciones declaradas 
%User--Inputs
N= 10; %Numero de subintervalos que dividen a [a,b]
a=0;
b=1;
f_x=@f2; %Declaracion de una funcion, Gilat Seccion 7.8       jjjjjjj
g_x=@f1; %Declaracion de una funcion, Gilat Seccion 7.8       jjjjjjj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
h=(b-a)/N;
x=linspace(a,b,N+1);%Crea un vector FILA de N puntos que van de [a,b]
                  % y cada punto esta separado por h=(b-a)/(N-1)
                  % Ver libro Gilat-Matlab pg 38        
y=f_x(x) + g_x(x);
plot(x,y,'-r*') %Ver Gilat seccion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7


%Definicion de la funcion f1(x)
function [val] = f1(x)
    val = x.^2;         %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
end


%Definicion de la funcion f2(x)
function [val] = f2(x)
    val = sin(pi*x);         %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
end
