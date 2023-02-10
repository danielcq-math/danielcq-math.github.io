%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
%Input
nI= 10; %Numero de subintervalos que dividen a [a,b]
a=0;
b=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
h=(b-a)/nI;       %longitud de cada subintervalo
x=zeros(nI+1,1);  %Creamos un vector columna con nI+1 entradas. Ver Gilat seccion 2.2.1
f_x=zeros(nI+1,1);
x(1)=a;       %primer elemento. Los subindices empiezasn desde 1 en Matlab
f_x(1)=x(1)^2;%primer elemento
%Ciclo para construir arreglos
for i=1:nI  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
    x(i+1)=a + i*h;
    f_x(i+1)=x(i+1)^2;
    
end

plot(x,f_x,'-r*'); %Gilat seccion 5.1