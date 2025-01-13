import numpy as np

from ngsolve import *
import netgen.gui

from cross import create_cross_mesh
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k

ngsglobals.msg_level = 1


nel = 12
#it is certainly smaller than this
n_glass = 1.55
min_wavelength = 1

w_domain = 0.8
h_domain = 0.7
w_cross = 0.05
l_cross = 0.5
h_silver = 0.3
h_sub = 0.1
h_pml = 0.1
el_ag = min_wavelength/nel
el_air = min_wavelength/nel
el_sub = n_glass*min_wavelength/nel
d_cross = 0.2

mesh = create_cross_mesh(w_domain,h_domain,w_cross,l_cross,d_cross,
                         h_silver,h_sub,h_pml,el_ag,el_air,el_sub)

#name of all 3d domains
vollist = {"air":1, "ag":2, "sub":3, "pml_air":4,
           "pml_sub":5}

#name of all 2d domains
surflist = {
    'air_port': 1,
    'sub_port' : 2,
    'dirichlet' : 3,
    }
count = len(surflist)+1
for name in vollist.keys():
    surflist[name+"_x_min"] = count
    count+=1
    surflist[name+"_x_max"] = count
    count+=1
    surflist[name+"_y_min"] = count
    count+=1
    surflist[name+"_y_max"] = count
    count+=1

print(surflist)
# checking surface domains

surf = mesh.BoundaryCF(surflist, -1)
gu = GridFunction(H1(mesh), name='surfs')
gu.Set(surf, definedon=~mesh.Boundaries(''))
Draw(gu)
print("Checking 2D domains... press Enter to continue")
input()

#checking volume domains

vol = CoefficientFunction([vollist.get(m,0) for m in mesh.GetMaterials()])
Draw(vol,mesh,"vols",draw_surf=False)
print("Checking 3D domains... press Enter to continue")
input()


wl_list = np.arange(1,1.1,0.005)
porder = 3

for wl in wl_list:
    ag_coeff = (ag_n(wl)-ag_k(wl)*1j)**2
    #for now, silver instead of glass
    sub_coeff = (ag_n(wl)-ag_k(wl)*1j)**2
    air_coeff = 1

    # defining permittivity for volume domains
    erlist = {
        'pml_air' : 1,
        'air' : 1,
        'ag' : (ag_n(wl)-ag_k(wl)*1j)**2,
        #for now, substrate is silver
        'sub' : (ag_n(wl)-ag_k(wl)*1j)**2,
        'pml_sub' : (ag_n(wl)-ag_k(wl)*1j)**2,
        }

    er = CoefficientFunction([erlist.get(m, 0) for m in mesh.GetMaterials()])
    #non-magnetic materials only
    ur = CoefficientFunction(1)

    #now we define the PML regions
    alphapml = 2.0j/h_pml

    boxmin = [-w_domain, -w_domain, 0]
    boxmax = [w_domain, w_domain, h_domain]
    mesh.SetPML(pml.Cartesian(mins=boxmin, maxs=boxmax, alpha=alphapml),"pml_air|pml_sub")

    fes = Periodic(HCurl(mesh, order=porder, complex=True,dirichlet="dirichlet"))

    u, v = fes.TnT()

    kzero = 2*pi/wl
    a = BilinearForm(fes,hermitian=False)
    a += ((1./ur) * curl(u) * curl(v) - kzero**2 * er * u * v)*dx

    f = LinearForm(fes)

    #our system is excited at the subdomain air port by a plane wave in the x direction
    #which is a constant vector field propagating with the propagation constant k0
    e_in = CoefficientFunction((0, 1))
    print(e_in)
    f += (2j*kzero * e_in * v.Trace()).ds("air_port")
    print("assembling system with ndofs {}".format(sum(fes.FreeDofs())))

    with TaskManager():
        a.Assemble()
        f.Assemble()

    # the solution field 
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
    


    # plot the solution (netgen-gui only)
    Draw (gfu)

    input()
