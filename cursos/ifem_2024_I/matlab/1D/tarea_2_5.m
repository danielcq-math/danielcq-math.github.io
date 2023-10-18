%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Teoría y Práctica Elementos Finitos
% Lic. en Matematicas-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descripcion: Dada la funcion u. Resolvemos el problema eliptico:
%                -d^2u/dx^2 = f en (0,1) con condiciones Dirichlet
%              utilizando
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
f_rhs=@f_exact; % Declaracion de una funcion handler (Gilat Ch 7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construimos un vector que obtiene el numero total de puntos en el mallado
%para cada ciclo de refinamiento 
nI_approx_vec=zeros(1,n_global_cicles);%contiene el numero inicial de nI_approx y su refinamiento
nI_approx_vec(1)=nI_approx_init; %El primer ciclo inicial (no se refina)

for k=2:n_global_cicles
  nI_approx_vec(k)=2*nI_approx_vec(k-1); % El siguiente contiene el doble de puntos que el anterior
end
 

L2_error_vec=zeros(1,n_global_cicles); % Vector que contiene el error en L2 para cada ciclo
% ciclos de refinamiento
for k=1:n_global_cicles

    nI_approx=nI_approx_vec(k);        % numero de subintervalos  en el mallado para aproximar
    % Vector que guarda los grados de libertdad de u(x), i.e, los valores
   
    nodes=linspace(a,b,nI_approx+1); %nodes de la malla, vector fila
    n_dofs=length(nodes); % numero de DOFs (grados de libertad, numero de nodos)
    A_global=zeros(n_dofs,n_dofs); %Matriz de difusion
    b_global=zeros(n_dofs,1);% vector RHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Proceso de ensamble
    %Ciclo para cada subintervalo
     for i=1:nI_approx
       %Operaciones locales
       b_local=int_phi_f_rhs(nodes,i); % b(1) tiene la contribucion para phi_i, 
                                       % b(2) tiene la contribucion para phi_(i+1)
       % Integrales locales de:  Dphi{i,2}*Dphi{i+1,1} en [i,i+1]
       % como la particion es uniforme podemos calcular esta integrales analiticamente:
       h=1/nI_approx;% h=(x_(i+1) - x_i))
       A_local(1,1)=1/h; %integral de Dphi{i,2}*Dphi{i,2}
       A_local(2,2)=1/h; %integral de Dphi{i+1,1}*Dphi{i+1,1}
       A_local(1,2)=-1/h;%integral de Dphi{i,2}*Dphi{i+1,1}
       A_local(2,1)=-1/h;%integral de Dphi{i+1,1}*Dphi{i,2}


       %Ensamble de b
       b_global(i)  = b_global(i)+ b_local(1);
       b_global(i+1)= b_global(i+1)+ b_local(2);

       %Ensamble de A
       A_global(i,i)    =A_global(i,i)+A_local(1,1);
       A_global(i+1,i+1)=A_global(i+1,i+1)+A_local(2,2);
       A_global(i,i+1)  =A_global(i,i+1)+A_local(1,2);
       A_global(i+1,i)  =A_global(i+1,i)+A_local(2,1);
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Manejo de las condiciones de Dirichlet
    %Condicion en x=a
    b_global(1)=u_exact_f(nodes(1));
    A_global(1,:)=zeros(1,n_dofs);
    A_global(1,1)=1;
    %Condicion en x=b
    b_global(n_dofs)=u_exact_f(nodes(n_dofs));
    A_global(n_dofs,:)=zeros(1,n_dofs);
    A_global(n_dofs,n_dofs)=1;

    %Resolver sistema
    uh_dofs=linsolve(A_global,b_global); 
    %uh_dofs=u_exact_f(nodes)'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculo del eror u - u_h en la norma L2
    L2_error_f=@L2_error;
    error=0;
    %Ciclo para cada subintervalo
    for i=1:nI_approx  %Ver Gilat-Matlab seccion 6.4 para ver ciclos for
        error= error+ L2_error_f(nodes, uh_dofs,i);
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
output_table=[nI_approx_vec' L2_error_vec' err_rate'] ; % Crea tabla con vectores columna
% %Ver Gilat-Matlab seccion 4.3 para el comando <disp>
disp("Tabla:"); %Imprime en la terminal                                                 
disp("nI_approx'  L2_enorm'  L2_erate'"); %Imprime en la terminal 
disp(output_table);%Imprime en la terminal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funciones definidas por el usuario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
%   Ver Libro de Gilat-Matlab Capitulo 7
% ** Esta funcion la define el usuario**

%Definicion (usuario) de la funcion  u_exact
function [val] = u_exact(x)
    val = sin(4*pi*x); 
end

%Definicion (usuario) de la funcion  f (RHS)
function [val] = f_exact(x)
    val = (4*pi).^2*sin(4*pi*x); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Funciones definidas por el programador
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: eval_interp
% Input: 
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i      : el índice del intervalo [nodes(i), nodes(i+1)]
%              
% Descripcion:
%               Evalua en el punto x el interpolador en el intervalo
%               [nodes(i),nodes(i+1)],
%               el cual
%               define como la combicion lineal del par de funciones con índices 
%               i e i+1 utilizando a su vez   los grados de libertad
%               dofs.
%               Precaucion: Para hacer esta función mucho más compacta
%               se asume que x efectivamente esta en el
%               intervalo [nodes(i), nodes(i+1)]
function [val] = eval_interp(x,nodes, dofs,i)

        eps=1e-12; %Tolerancia
        %Checamos si x esta en el intervalo donde hat_i(x) es cero
        if((x+eps<nodes(i)) || (x-eps > nodes(i+1))) % Ver Gilat sec. 6.1 & 6.2
           disp("Error: x is not defined in eval_hat_functions_pair ") ;%Gilat sec 4.3.1
           assert(false); %Stops program
        end
 
       val= dofs(i)*(1- (x- nodes(i))/(nodes(i+1)-nodes(i)));% phi_{i,2}
       val= val + dofs(i+1)*(x-nodes(i))/(nodes(i+1)-nodes(i));% phi_{i+1,1}
    
end
% Function: int_phi_f_rhs
% Input: 
%       x:      punto a evaluar
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes
%       i         : el índice del intervalo [nodes(i), nodes(i+1)] 
% Output:
%      val:  val(1) tiene la contribucion de integral de phi_i_2*f_rhs
%            val(2) tiene la contribucion de integral de phi_(i+1)_1*f_rhs
%              
% Descripcion:
%               Aproxima las integrales de phi_hat*f_rhs 
%               en el intervalo
%               [nodes(i),nodes(i+1)] utilizando la regla de Simpson
%                 

function [val] = int_phi_f_rhs(nodes,i)
    %Regla de Simpson
    %Puntos extemos de la integral son x0,x2    
            x0 =nodes(i);
            x2 =nodes(i+1);
            k=0.5*(x2-x0);
            x1 = x0+k;
            gpts=[x0 x1 x2];
            gws=k/3.0*[1.0 4.0 1.0];
            val=zeros(2,1);
            for l=1:3
               x=gpts(l);
               w=gws(l);
               phi_i_2=1-(x-nodes(i))/(nodes(i+1)-nodes(i));% phi_{i,2}
               phi_ip1_1=(x-nodes(i))/(nodes(i+1)-nodes(i));% phi_{i+1,1}
               val(1)=val(1) + w*f_exact(x)*phi_i_2;
               val(2)=val(2) + w*f_exact(x)*phi_ip1_1;
            end

end

% Function: L2-error
%       nodes:  los nodos de la particion [a,b]
%       dofs:   los grados de libertad en nodes       
%       i      : el índice de la hat_función  de la izquierda
% Descripcion:
%       Evalua la integral (u-uI)^2 en el intervalo  (nodes(i),nodes(i+1))
%       utilizando  la regla de Simpson. Ver Libro de Burden-Faires pg.196 
function [val] = L2_error(nodes, dofs,i)


            %Puntos extemos de la integral son x0,x2    
            x0 =nodes(i);
            x2 =nodes(i+1);
            h=0.5*(x2-x0);
            x1 = x0+h;
          
            val= (u_exact(x0) - eval_interp(x0,nodes, dofs,i))^2;
            val= val + 4*(u_exact(x1) - eval_interp(x1,nodes, dofs,i))^2;
            val= val + (u_exact(x2) - eval_interp(x2,nodes, dofs,i))^2;
            val = h/3.*val;
end

