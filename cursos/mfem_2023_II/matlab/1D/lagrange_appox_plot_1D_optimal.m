%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Posgrado en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada una función u(x)
%              obtenemos su interpolador u_h(x) utilizando
%              las funciones hat_j(x) lineales por pedazos, i.e., 
%               funciones de lagrange P1.
%Ouput: gráfica de u(x) y u_h(x) (programacion de forma optima**)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx=10;        % numero de subintervalos  en el mallado para aproximar
u_exact_f=@u_exact; % Declaracion de una funcion handler (Gilat Ch 7)
nI_plot= 2*nI_approx;       % numero total de subintervalos en el mallado para graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vector que guarda los grados de libertdad de u(x), i.e, los valores
% u(xi)
nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
n_dofs=length(nodes); % numero de DOFs (degrees of freedom)
dofs=u_exact_f(nodes)';  %vector columna (transpuesto)
hat_f_pair=@hat_functions_pair;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables para Graficar 
h_plot=(b-a)/nI_plot;
x=zeros(nI_plot+1,1);       %Dominio de grafica
uh=zeros(nI_plot+1,1);       %vector para plotear en el rango de la grafica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grafica de u_h, ** forma  optima **
x(1)=a;         %primer elemento  ** los indices empiezan desde 1 

% Recordar que tenemos hat_functions lineales. Entonces para evaluar
% en un punto x  el span de hat_functions es suficiente
% encontrar las dos hat_functions que son distintas de cero en x.
% Por ejemplo para x en el intervalo semi-abierto [nodes(1), nodes(2)),
% las hat_functions con indices 1  y 2 son
% las únicas que son distintas de cero. 
% Por lo tanto definimos los 
% los indices de la hat_functions para este fin:
nzero_hat_idxs=zeros(2,1);
nzero_hat_idxs(1)=1;
nzero_hat_idxs(2)=2;

uh(1)=hat_f_pair(a,nodes,dofs,1);  %primer elemento
%Ciclo para construir arreglos
for i=1:nI_plot  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
    x(i+1)=a + i*h_plot;
         % x debe estar en el intervalo
         %   [   nodes(nzero_hat_idxs(1)),nodes(nzero_hat_idxs(2)) )
        uh(i+1)= hat_f_pair(x(i+1),nodes,dofs,nzero_hat_idxs(1));

    % Checar si tenemos que cambiar los nzero_hat_idxs(),i.e,
    %  x ya no en el intervalo [  nodes(nzero_hat_idxs(1)),nodes(nzero_hat_idxs(2)) )
    if(x(i+1)>=nodes(nzero_hat_idxs(2))) 
        nzero_hat_idxs(1)= nzero_hat_idxs(2);
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

% Hat function pair
% Input:
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i      : el índice de la hat_función  decla izquierda
%              
% Descripcion:
%               Evalua en el punto x el par de funciones con índices 
%               i e i+1 utilizando a su vez los  los grados de libertad
%               dofs.
%               Precaucion: Para hacer esta función mucho más compacta
%               se asume que x efectivamente esta en el
%               intervalo [nodes(i), nodes(i+1)]
function [val] = hat_functions_pair(x,nodes, dofs,i)

        eps=1e-6; %Tolerancia
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x+eps<nodes(i)) || (x-eps > nodes(i+1))) % Ver Gilat sec. 6.1 & 6.2
           disp("Error: x is not defined in hat_functions_pair ") ; %Gilat sec 4.3.1
           assert(false); %Stops program
        end
 
       val= dofs(i)*(1- (x- nodes(i))/(nodes(i+1)-nodes(i)));
       val= val + dofs(i+1)*(x-nodes(i))/(nodes(i+1)-nodes(i));
    
end