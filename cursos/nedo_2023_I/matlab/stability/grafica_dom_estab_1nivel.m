%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Facultad de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ** Grafica de dominios de estabilidad para métodos numéricos de un paso/nivel **
%
%
% Representamos un numero complejo z como z = x +i*y donde x,y son numeros
% reales, es decir, la variable x es la parte real y la variable y es la parte imaginaria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Parametros del usuario
%Seleccionar uno definicion de Dom (ver definiciones al final del archivo)
%Dom=@DomEulerEx_f;    %Metodo Euler explicito
%Dom=@DomEulerImp_f;   %Metodo de Euler implicito
Dom=@DomTaylorEx2_f;   %Metodo  de Taylor explicito de orden r=2
%Dom=@DomTrapecio_f;      %Metodo  del Trapecio
%Dominio para graficar
x_a=-2; x_b=2;
y_a=-2; y_b=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fronteraDom = @(x,y) 1 - Dom(x,y); %  Frontera del Dominio de estabilidad: {{Dom} = 1} 
fimplicit(fronteraDom,[x_a x_b y_a y_b],'r','LineWidth',2);

 hold on
 dentroDom = @(x,y) 0.5 - Dom(x,y); %  Algunos Puntos dentro del Dominio de estabilidad:  {{ Dom} = 0.5}
 fimplicit(dentroDom,[x_a x_b y_a y_b],'--g');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **La definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7

%Definicion de la funcion Dom para el Met de Euler Expplicito
function [val] = DomEulerEx_f(x,y)
   %   | 1 + z |,  DominioDeEstabilidad = {{ Dom} < 1} 
   %   z = x + i*y
       val =  sqrt( (1 + x).^2 + y.^2); 
end


%Definicion de la funcion Dom para el Met de Euler Implicito
function [val] = DomEulerImp_f(x,y)
    % | 1/(1 - z) |   DominioDeEstabilidad = {{ Dom} < 1}
    %   z = x + i*y
    val =   1/sqrt( (1 - x).^2 + y.^2); 
end


%Definicion de la funcion Dom para el M. Explicito de Taylor de orden r =2
function [val] = DomTaylorEx2_f(x,y)
    % | 1 + z + z^2 |   DominioDeEstabilidad = {{ Dom} < 1}
    %   z = x + i*y
    val =   sqrt((1+ x + 0.5*x.^2 - 0.5*y.^2).^2 + (y+x.*y).^2); 
end


%Definicion de la funcion Dom para el M.  del Trapecio
function [val] = DomTrapecio_f(x,y)
       %   | (1 + 0.5*z)/(1 - 0.5*z)   DominioDeEstabilidad = {{ Dom} < 1}
       %   z = x + i*y
    val =  sqrt((1-0.25*x.^2-0.25*y.^2)^2 + y.^2)/( (1 - 0.5*x).^2 + 0.25*y.^2); 
end


