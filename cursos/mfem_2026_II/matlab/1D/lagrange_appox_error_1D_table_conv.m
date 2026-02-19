%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Posgrado en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada una función u(x)
%              obtenemos su interpolador u_h(x) utilizando
%              las funciones hat_j(x) lineales por pedazos, i.e., 
%              funciones de lagrange P1.
%              El calculamos el error (u - uh) en la norma L2
%              para diferentes refinaciones del intervalo
%Ouput:    Tabla de convergencia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametros de Usuario
%interval endpoints of [a,b]
a=0;
b=1;
nI_approx_init=10;      % numero inicial del total de sub-intervalos en el mallado
n_global_cicles=5;             % numero de ciclos que corre el algoritmo numérico
u_exact_f=@u_exact; % Declaracion de una funcion handler (Gilat Ch 7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construimos un vector que obtiene el numero total de puntos en el mallado
%para cada ciclo de refinamiento 
nI_approx_vec=zeros(1,n_global_cicles);
nI_approx_vec(1)=nI_approx_init; %El primer ciclo no se refina

for k=2:n_global_cicles
  nI_approx_vec(k)=2*nI_approx_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end
 
L2_error_vec=zeros(1,n_global_cicles); % Vector que contiene el error en L2 para cada ciclo
% ciclos de refinamiento
for k=1:n_global_cicles

    nI_approx=nI_approx_vec(k);        % numero de subintervalos  en el mallado para aproximar
    % Vector que guarda los grados de libertdad de u(x), i.e, los valores
    % u(xi)
    nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
    n_dofs=length(nodes); % numero de DOFs (degrees of freedom)
    dofs=u_exact_f(nodes)';  %vector columna (transpuesto)
    L2_error_f=@L2_error;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma L2
    error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
        error= error+ L2_error_f(nodes, dofs,i);
    end
    error= sqrt(error);
    L2_error_vec(k)=error; %Guardar el valor del error L2
    output_text=['El error L2 para nI_approx=', num2str(nI_approx), ' de subintervalos es: ', num2str(L2_error_vec(k))];
    disp(output_text);
end %acaba el ciclo de refinamientos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% Crear tabla con valores: N| max(y-val - y_hal)| tasa_de_decrecimiento
disp("##### Final del algoritmo #########"); %Imprime en la terminal 
% decrecimiento del error:= log(err_max(i+1)/err_max(i))*log(2)
err_rate=zeros(1,n_global_cicles); %tasa de decrecimiento del error
err_rate(1)=1;
for i=2:n_global_cicles
    err_rate(i)=log(L2_error_vec(i)/L2_error_vec(i-1))/log(1/2);
end
output_table=[nI_approx_vec' err_rate'] ; % Crea tabla con vectores columna
% %Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla: nI_approx'  L2_err_rate'"); %Imprime en la terminal                                                 
disp(output_table);%Imprime en la terminal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
%   Ver Libro de Gilat-Matlab Capitulo 7
% ** Esta funcion la define el usuario**

%Definicion (usuario) de la funcion  u_exact
function [val] = u_exact(x)
    val = cos(4*pi*x); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functio: eval_hat_function_pair
% Input: 
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i      : el índice de la hat_función  de la izquierda
%              
% Descripcion:
%               Evalua en el punto x el par de funciones con índices 
%               i e i+1 utilizando a su vez los  los grados de libertad
%               dofs.
%               Precaucion: Para hacer esta función mucho más compacta
%               se asume que x efectivamente esta en el
%               intervalo [nodes(i), nodes(i+1)]
function [val] = eval_hat_function_pair(x,nodes, dofs,i)

        eps=1e-6; %Tolerancia
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x+eps<nodes(i)) || (x-eps > nodes(i+1))) % Ver Gilat sec. 6.1 & 6.2
           disp("Error: x is not defined in hat_functions_pair ") ;%Gilat sec 4.3.1
           assert(false); %Stops program
        end
 
       val= dofs(i)*(1- (x- nodes(i))/(nodes(i+1)-nodes(i)));
       val= val + dofs(i+1)*(x-nodes(i))/(nodes(i+1)-nodes(i));
    
end

% Function: L2-error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  del a izquierda
% Descripcion:
%       Evalua la integral (u-uh)^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  la regla del Simpson. Ver Libro de Burden-Faires pg.196 
function [val] = L2_error(nodes, dofs,i)

            x0 =nodes(i);
            x2 =nodes(i+1);
            h=0.5*(x2-x0);
            x1 = x0+h;
          
            val= (u_exact(x0) - eval_hat_function_pair(x0,nodes, dofs,i))^2;
            val= val + 4*(u_exact(x1) - eval_hat_function_pair(x1,nodes, dofs,i))^2;
            val= val + (u_exact(x2) - eval_hat_function_pair(x2,nodes, dofs,i))^2;
            val = h/3.*val;
end
