%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
% utilizando aritmetica vectorial de matlab
%Input
nI= 10;  %Numero de subintervalos que dividen a [a,b]
a=0;
b=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
h=(b-a)/nI;        %longitud de cada subintervalo
x=linspace(a,b,nI+1);%Crea un vector FILA de nI + 1 puntos que van de [a,b]
                  % y cada punto esta separado por h=(b-a)/(N-1)
                  % Ver libro Gilat-Matlab pg 38 
f_x=x.^2; %Aritmetica vectorial, libro Gilat-Matlab seccion 3.4
plot(x,f_x,'-r*')