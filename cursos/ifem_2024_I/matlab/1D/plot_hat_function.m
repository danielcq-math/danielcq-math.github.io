%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Curso Introduccion a los Elementos Finitos
% Lic. en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Descripcion:  Grafica 
%               la funcion hat_i(x)  por pedazos, i.e., 
%               funciones de lagrange de primer orden.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx=10;        % numero de subintervalos  en el mallado para definir las funciones hat
index_i= 8;         %indice de la funcion hat_i(x) a graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
n_dofs=length(nodes);    % numero de nodos (degrees of freedom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros para Graficar hat_i(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%El dominio y el numero de puntos a graficar coincide con el numero de nodos
n_plot=n_dofs; 
x_dominio=linspace(nodes(1),nodes(end),n_plot);%Dominio para graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grafica solo el dominio de hat_i(x) donde es distinta de cero
%x_dominio=linspace(nodes(index_i-1),nodes(index_i+1),n_plot);%Dominio para graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graficar hat_i(x)

y=zeros(n_plot,1); %vector rango
for i=1:n_plot
    y(i)=eval_hat_function(x_dominio(i),nodes,index_i); %Vector rango
end
plot(x_dominio,y,'-r*') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
%   Ver Libro de Gilat-Matlab Capitulo 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funciones definidas por el programador
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: eval_hat_function
% Input: 
%       x:      punto a evaluar, 
%%              x tiene que ser un numero
%       nodes:  los nodos de la particion [a,b]
%       i    :  el índice de la hat_función que se quiere evaluar
%              
% Descripcion:
%               Evalua en el punto x la  funcion hat con índice i.
%             

function [val] = eval_hat_function(x,nodes,i)

        val=0; %valor por default

        eps=1e-12; %Tolerancia
      
 
       if((x+eps> nodes(i-1)) && (x-eps < nodes(i))) 
           val = (x-nodes(i-1))/(nodes(i)-nodes(i-1)); %phi_{i,1}
       elseif((x+eps> nodes(i)) && (x-eps < nodes(i+1))) 
           val= (1- (x- nodes(i))/(nodes(i+1)-nodes(i)));%phi_{i,2}
       end
    
end


