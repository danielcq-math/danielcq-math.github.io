%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Facultad de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grafica de dominios de estabilidad
% Representamos un numero complejo z como z = x +i*y donde x,y son numeros
% reales, es decir, la variable x es la parte real y la variable y es la parte imaginaria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Euler explicito
DomEuler= @(x,y) sqrt( (1 + x).^2 + y.^2); % | 1 + z |,  DominioEuler = {{ DomEuler} < 1}

fronteraDom = @(x,y) 1 - DomEuler(x,y); %  {{ DomEuler} = 1}
fimplicit(fronteraDom,[-2 2 -1 1],'r','LineWidth',2);

 hold on
 dentroDom = @(x,y) 0.5 - DomEuler(x,y); % {{ DomEuler} = -5}
 fimplicit(dentroDom,[-2 2 -1 1],'--g');
 hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%