######################################################
##% Curso de Teoría y Práctica Elementos Finitos
#% Posgrado en Matematicas-UNAM-CdMx
#% Prof. Daniel Castañon Quiroz. daniel.castanon@iimas.unam.mx
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Descripcion: Dada una función u(x)
#%              obtenemos su interpolador u_h(x) utilizando
#%              las funciones hat_j(x) lineales por pedazos, i.e., 
#%              funciones de lagrange P1.
#%              El calculamos el error (u - uh) en la norma L2
#%              para diferentes refinaciones del intervalo
#%Ouput:     Error en H1 y L2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% Parametros de Usuario
#%% hmax, orden del EF de Lagrange
####################################################
#This example is taken from ngsolve-itutorials and the iFEM/primal course


from ngsolve import *
#from ngsolve.webgui import Draw

################################################
#Parametros de usuario
hmax=0.05
FE_order=1  #Grado del EF de Lagrange
#################################################

#unit square [0,1]^2
mesh = Mesh(unit_square.GenerateMesh(maxh=hmax))
#Draw (mesh);
fes = H1(mesh, order=FE_order, dirichlet="left|right|bottom|top")
print ("ndof =", fes.ndof)
print (mesh.GetBoundaries())
u = fes.TrialFunction()
v = fes.TestFunction()

f = LinearForm(fes)
f += 32*(y*(1-y)+x*(1-x))*v*dx

a = BilinearForm(fes)
a += grad(u)*grad(v)*dx   

a.Assemble()
f.Assemble();
#print(f.vec)
#print(a.mat)
######################################################
#Resuelve el sistema lineal
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
#######################################################


exact = 16*x*(1-x)*y*(1-y)
D_exact = (16*y*(1-y)*(1-2*x),16*x*(1-x)*(1-2*y))

#Print results

print ("h_max:", hmax)
print ("FE_order:", FE_order)
print ("L2-error:", sqrt(Integrate((gfu-exact)**2, mesh)))
print ("H1-error:", sqrt(Integrate((grad(gfu)-D_exact)**2, mesh)))

