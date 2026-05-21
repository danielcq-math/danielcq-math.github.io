####################################################################################################
##% Curso de Teoría y Práctica Elementos Finitos
#% Posgrado en Matematicas-UNAM-CdMx
#% Prof. Daniel Castañon Quiroz. danielcq@iimas.unam.mx
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#################################
#% Descripcion: Resuelve el problema de Stokes con el par estable Taylor-Hood P2-P1. 
#  Condiciones  Dirichlet nulas
#%Ouput:     Error en H1 y L2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###################################################
#################################################################################################
#This example is inspired  from ngsolve-itutorials Section 2.6
# and the iFEM/primal course Ch 32
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###################################################
#%% Parametros de Usuario
hmax=0.05
mu=1 #viscosidad
#########################################################################################
#
from ngsolve import *
#from ngsolve.webgui import Draw
from math import pi
from netgen.geom2d import SplineGeometry

#unit square [0,1]^2
geo = SplineGeometry()
geo.AddRectangle((0,0),(1,1),bc="rectangle")

###########################################################
# exact velocity
ex_vel =  ( 0.5*sin(2*pi*x)*cos(2*pi*y),
           -0.5*cos(2*pi*x)*sin(2*pi*y))
# exact pressure
ex_pre = x -0.5
# source term
src =  ( 4*mu*pi**2*sin(2*pi*x)*cos(2*pi*y) + 1,
         -4*mu*pi**2*cos(2*pi*x)*sin(2*pi*y))


##########################################################

mesh = Mesh(geo.GenerateMesh(maxh=hmax))
_order_=1
mesh.Curve(order=_order_)
print("ElementBoundary=",mesh.GetBoundaries())
print("NumEles=",mesh.ne)
print("NumEdges=",mesh.nedge)
print("NumVertex=",mesh.nv)
 
#Taylor Hood Pair
V = VectorH1(mesh, order=2, dirichlet="rectangle")
Q = H1(mesh, order=1)

X = V*Q

u_exact = CoefficientFunction(ex_vel)

def SolveStokes(X):
    u,p = X.TrialFunction()
    v,q = X.TestFunction()
    a = BilinearForm(X)
    a += (InnerProduct(grad(u),grad(v))-div(u)*q-div(v)*p)*dx
    a += -1e-10*p*q*dx #Necesary in practice
    a.Assemble()
    # a.mat[X.components[0].ndof,X.components[0].ndof]+=-1  # fixing pressure in one point
    
    
    f = LinearForm(X)
    source = CoefficientFunction(src)
    f += InnerProduct(source,v)*dx
    f.Assemble()
    
    
    #Exact Dirichelt BC for velocity
    gfu = GridFunction(X)
    gfu.components[0].Set(u_exact, definedon=mesh.Boundaries("rectangle"))    
    
    
    ## solve the linear system
    ## Check how to set inhomogeneus Dirichlet Conditions
    res = f.vec.CreateVector()
    res.data = f.vec - a.mat*gfu.vec
    inv_stokes = a.mat.Inverse(X.FreeDofs())
    gfu.vec.data += inv_stokes * res
    #Draw(gfu)

    return gfu


gfu=SolveStokes(X)

#uh, ph = gfu.components
uh = CoefficientFunction(gfu.components[0])
ph = CoefficientFunction(gfu.components[1])
p_exact = CoefficientFunction(ex_pre)

err_L2_u = sqrt(Integrate(InnerProduct(uh - u_exact, uh - u_exact), mesh))
err_p = sqrt(Integrate((p_exact-ph)**2, mesh))

print ("h_max:", hmax)
print ("L2-vel error:", err_L2_u)
print ("L2-presure error:", err_p)

