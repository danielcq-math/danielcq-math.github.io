%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Posgrado en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada una función u(x)
%              obtenemos su interpolador u_h(x) utilizando
%              las funciones hat_j(x) lineales por pedazos, i.e., 
%               funciones de lagrange P1.
%Ouput: gráfica de u(x) y u_h(x) (programacion de forma no optima**)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx=10;        % numero de subintervalos  en el mallado para aproximar
u_exact_f=@u_exact; % Declaracion de una funcion handler (Gilat Ch 7)
nI_plot= 4*nI_approx;       % numero total de subintervalos en el mallado para graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vector que guarda los grados de libertdad de u(x), i.e, los valores
% u(xi)
nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
n_dofs=length(nodes); % numero de DOFs (degrees of freedom)
dofs=u_exact_f(nodes)';  %vector columna (transpuesto)
hat_f=@hat_function;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables para Graficar 
h_plot=(b-a)/nI_plot;
x=zeros(nI_plot+1,1);       %Dominio de grafica
y=zeros(nI_plot+1,1);       %Rango
uh=zeros(nI_plot+1,1);       %vector para plotear en el rango de la grafica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %Grafica una hat function
%% Se puede comentar todo el bloque si no se quiere graficar la funcion hat

% %% elegir un indice de hat_function 
% hat_idx=n_dofs; %Probar desde 1 hast n_ndofs
% x(1)=a;         %primer elemento  ** los indices empiezan desde 1 
% y(1)=hat_f(a,nodes,hat_idx);  %primer elemento
% %Ciclo para construir arreglos
% for i=1:nI_plot  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
%     x(i+1)=a + i*h_plot;
%     y(i+1)=hat_f(x(i+1),nodes,hat_idx);    
% end
% %Grafica la hat function con indice hat_idx
%  plot(x,y,'-r'); %Ver Gilat seccion 5.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grafica de u_h, ** forma no optima **
x(1)=a;         %primer elemento  ** los indices empiezan desde 1 
uh(1)=hat_f(a,nodes,1);  %primer elemento
%Ciclo para construir arreglos
for i=1:nI_plot  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
    x(i+1)=a + i*h_plot;
    uh(i+1)=0;
    for j=1:n_dofs %loop de todas la hat functions **no optimo)
        uh(i+1)= uh(i+1) + dofs(j)*hat_f(x(i+1),nodes,j);  
    end
end
u=u_exact_f(x); %valor de la funcion smooth para comparar
%Grafica uh (rojo), u (azul) y los valores nodales con *
plot(x,uh,'-r',x,u,'b',nodes,dofs,'*r'); %Ver Gilat seccion 5.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
%   Ver Libro de Gilat-Matlab Capitulo 7
% ** Esta funcion la define el usuario**

%Definicion (usuario) de la funcion  u_exact
function [val] = u_exact(x)
    val = cos(4*pi*x); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function hat
%Input: x:      donde se evalua la funcion
%       nodes:  los nodos de la particion [a,b]
%       i    :  el índice de la función que se quiere evaluar
function [val] = hat_function(x,nodes,i)


        val=0; %valor por default

        n_nodes=length(nodes);
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x<nodes(1)) || (x > nodes(n_nodes))) % Ver Gilat sec. 6.1 & 6.2
           return; 
        end

        %Para i > 2
        if((i>1) && (x < nodes(i-1)))
            return; 
        end

        %Para hat_(n_nodes + 1)
        if( (i<n_nodes) && (x > nodes(i+1))) 
            return; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % x esta en el dominio de hat_i donde es distinta de cero
        
        if(i==1) %first hat function
            val= 1-(x-nodes(1))/(nodes(2)-nodes(1));
            return; 
        end


        if(i==n_nodes) %last hat function
            val= (x-nodes(n_nodes-1))/(nodes(n_nodes)-nodes(n_nodes-1));
            return; 
        end


        %Other hat functions,i.e, 2 <= i <= n_nodes-1
        if((x > nodes(i-1)) && (x < nodes(i)))
             val= (x- nodes(i-1))/(nodes(i)-nodes(i-1));
             return; 
        end

        if((x >= nodes(i)) && (x < nodes(i+1))) %Need to put >=
             val= 1- (x- nodes(i))/(nodes(i+1)-nodes(i));
             return; 
        end
end