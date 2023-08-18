%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curso Intro FEM
% Faculta de de Ciencias-UNAM-CdMx
% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Indices y vectores
% Secciones 2.1 y 2.5 libro de Gilat-Matlab

%Creando un vector directamente escribiendo sus valores
x=[0 1.3 5];

%Usando el comando linspace
y=linspace(0,1,10);

%Usando el comando ones() y zeros()
z=ones(3,1);
w=zeros(1,5);

%Operaciones
x_trans= x'; %transpuesto
y2= x+ 3*x;


%Acceso
x(1); %primer elemento
x(2);  %segundo elemento
x(1:3); %Elementos del primer al tercero
x(1:2); %elementos del primer al segundo