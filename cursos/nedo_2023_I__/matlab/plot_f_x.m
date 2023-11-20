%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica de una funcion f(x) en el intervalo [a,b]
%Input
N= 10; %Numero de puntos
a=0;
b=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Init
h=(b-a)/N;
x=zeros(N,1);
f_x=zeros(N,1);
x(1)=a;
f_x(1)=x(1)^2;%primer elemento
%Ciclo para construir arreglos
for i=1:N  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
    x(i+1)=a + i*h;
    f_x(i+1)=x(i+1)^2;
    
end

plot(x,f_x,'-r*');