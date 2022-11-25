%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso de Solución Numérica de Ecuaciones Diferenciales Ordinarias
% Facultad de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ** Grafica de dominios de estabilidad para métodos numéricos MULTI-paso/nivel **
%
%
% Representamos un numero complejo z como z = x +i*y donde x,y son numeros
% reales, es decir, la variable x es la parte real y la variable y es la parte imaginaria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs del usuario
clearvars %limpia los valores de variables definidas anteriormente 
N_theta= 75; %numero de sub-intervalos que dividen el intervalo global [0,2pi]
r0= 0.975;   %valor radial para graficar puntos dentro del dominio de estabilidad

%Seleccinar solo un rho_f y un sigma_f:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metodo de Adams-Bashfort (Exp) de  2 niveles
rho_f=@rho_AdamsB2;
sigma_f=@sigma_AdamsB2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Metodo de Adams-Bashfort (Exp) de  4 niveles
%rho_f=@rho_AdamsB4;

%sigma_f=@sigma_AdamsB4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metodo de Adams-Moulton (Imp) de  2 niveles
%rho_f=@rho_AdamsM2;
%sigma_f=@sigma_AdamsM2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BDF2 implicito
%rho_f=@rho_BDF2;
%sigma_f=@sigma_BDF2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_v=linspace(0,2*pi,N_theta);  %comando que divide [0,2*pi] con Ntheta puntos

%frontera del DomEstab  
psi_bd= cos(theta_v)+i*sin(theta_v); %construimos z= exp(i*theta)= cos(theta) +i*sin(theta)
z_bd=rho_f(psi_bd)./sigma_f(psi_bd); %frontera

%Dentro del DomEstab  
psi_in= r0*(cos(theta_v)+i*sin(theta_v));%construimos z= r0*exp(i*theta)
z_in=rho_f(psi_in)./sigma_f(psi_in); %dentro


x_bd=real(z_bd); %parte real 
y_bd=imag(z_bd); %parte imaginaria
x_in=real(z_in);
y_in=imag(z_in);
plot(x_bd,y_bd,'r',x_in,y_in,'g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **Las definiciones de las funciones deben ir al final del archivo script
% Ver Libro de Gilat-Matlab Capitulo 7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adams-Bashfort m=2 niveles
function [val] = rho_AdamsB2(psi)
    val = psi.^2 -psi; 
end

function [val] = sigma_AdamsB2(psi)
    [val] = 3/2*psi - 1/2; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adams-Bashfort m=4 niveles
function [val] = rho_AdamsB4(psi)
    val = psi.^4 -psi.^3; 
end

function [val] = sigma_AdamsB4(psi)
    [val] = 55/24*psi.^3 - 59/24*psi.^2 + 37/24*psi - 9/24; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adams-Moulton m=2 niveles
function [val] = rho_AdamsM2(psi)
    val = psi.^2 -psi; 
end

function [val] = sigma_AdamsM2(psi)
    [val] = 5/12*psi.^2 + 2/3*psi - 1/12; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BDF2 implicito
function [val] = rho_BDF2(psi)
    val = 1.5*psi.^2 -2*psi+0.5; 
end

function [val] = sigma_BDF2(psi)
    [val] = psi.^2; 
end