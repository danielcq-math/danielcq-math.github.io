####################################################################################################
##% Curso de Teoría y Práctica Elementos Finitos
#% Posgrado en Matematicas-UNAM-CdMx
#% Prof. Daniel Castañon Quiroz. danielcq@iimas.unam.mx
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#################################
#% Descripcion: Resuelve el problema mixto de Poisson con E.F de Raviar-Thomas
#%Ouput:     Error en L2 del flujo (gradiente) y L2 del potencial
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###################################################
#%% Parametros de Usuario
#%% hmax, orden del EF de Raviar-Thomas-Nedelec
#################################################################################################
#This example is taken from https://jschoeberl.github.io/iFEM/mixed/secondorder.html
#Mixed problem for Poisson Eq. with zero Dirichlet conditions


from ngsolve import *
#from ngsolve.webgui import Draw

####################################################
## Parametros del usuario ###
hmax=0.025
Hdiv_order=1 #orden del espacio de RT

source= 32*(y*(1-y)+x*(1-x))
u_exact = 16*x*(1-x)*y*(1-y)
D_uexact = (16*y*(1-y)*(1-2*x),16*x*(1-x)*(1-2*y))


########################################################
#unit square [0,1]^2
mesh = Mesh(unit_square.GenerateMesh(maxh=hmax))
#Draw (mesh);
#RT elements for simplicial element
RT_FE_flag= True #Use RT FE or BDM FE

Sigma = HDiv(mesh, order=Hdiv_order, RT=RT_FE_flag)#If RT false then it is BDM
L2_order=Hdiv_order-1 #Default
if RT_FE_flag:
 L2_order=Hdiv_order

M = L2(mesh, order=L2_order)
X = Sigma*M

sigma,u = X.TrialFunction()
tau,v = X.TestFunction()

#Formular como BNB
a = BilinearForm( (sigma*tau+div(sigma)*v+div(tau)*u)*dx).Assemble()
f = LinearForm(-source*v*dx).Assemble()

print ("ndof =", X.ndof)
#print (mesh.GetBoundaries())
#print(f.vec)
#print(a.mat)
gfu = GridFunction(X)
gfu.vec.data = a.mat.Inverse(X.FreeDofs()) * f.vec

flux_h=CoefficientFunction(gfu.components[0])
uh=CoefficientFunction(gfu.components[1])

#Print results

print ("h_max:", hmax)
print ("FE_order:", Hdiv_order)
print ("L2-Flux-error:", sqrt(Integrate((flux_h-D_uexact)**2, mesh)))
print ("L2-error:", sqrt(Integrate((uh-u_exact)**2, mesh)))
