%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
% utilizando aritmetica vectorial de matlab
%User--Inputs
N= 10; %Numero de subintervalos que dividen a [a,b]
a=0;
b=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
h=(b-a)/N;
x=linspace(a,b,N+1);%Crea un vector FILA de N puntos que van de [a,b]
                  % y cada punto esta separado por h=(b-a)/(N-1)
                  % Ver libro Gilat-Matlab pg 38 
y=x.^2;         %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
plot(x,y,'-r*') %Ver Gilat seccion 5.1