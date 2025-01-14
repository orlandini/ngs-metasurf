from netgen.occ import *
from netgen.meshing import IdentificationType
from ngsolve import Mesh


def create_cross_mesh(w_domain, h_domain, w_cross, l_cross, d_cross,
                      h_silver,h_sub,h_pml,el_ag, el_air, el_sub):

    """
    Creates a mesh representing the unit cell of a metasurface consisting
    of Ag layer in which a cross was carved over a substrate
    """
    #initial height of cross
    h_cross = h_sub + h_silver - d_cross

    box = Box(Pnt(-w_domain/2,-w_domain/2,0), Pnt(w_domain/2,w_domain/2,h_domain))
    sub = Box(Pnt(-w_domain/2,-w_domain/2,0), Pnt(w_domain/2,w_domain/2,h_sub))
    ag_full = Box(Pnt(-w_domain/2,-w_domain/2,h_sub), Pnt(w_domain/2,w_domain/2,h_silver+h_sub))
    cross_1 = Box(Pnt(-w_cross/2,-l_cross/2,h_cross), Pnt(w_cross/2,l_cross/2,h_silver+h_sub))
    cross_2 = Box(Pnt(-l_cross/2,-w_cross/2,h_cross), Pnt(l_cross/2,w_cross/2,h_silver+h_sub))

    #pml regions
    pml_air = Box(Pnt(-w_domain/2,-w_domain/2,h_domain), Pnt(w_domain/2,w_domain/2,h_domain+h_pml))
    pml_sub = Box(Pnt(-w_domain/2,-w_domain/2,-h_pml), Pnt(w_domain/2,w_domain/2,0))

    #boolean operations
    cross = cross_1 + cross_2
    ag = ag_full - cross
    air = box - sub - ag

    #now for element sizes
    air.mat("air").maxh = el_air
    sub.mat("sub").maxh = el_sub
    ag.mat("ag").maxh = el_ag

    pml_air.mat('pml_air').maxh = el_air
    pml_sub.mat('pml_sub').maxh = el_sub


    
    domain_dict = {
        pml_air : "pml_air",
        air : "air",
        ag : "ag",
        sub : "sub",
        pml_sub : "pml_sub"
    }
    domain_list = list(domain_dict.keys())
    #useful for checking later
    for domain in domain_list:
        domain.faces.name = 'default'

    #setting the ports
    sub.faces.Min(Z).name = 'sub_port'
    #this does not work
    air.faces.Max(Z).name = 'air_port'
    #but this does
    pml_air.faces.Min(Z).name = 'air_port'


    
    #setting dirichlet BCs at PMLs
    pml_air.faces.Max(Z).name = 'dirichlet'
    pml_sub.faces.Min(Z).name = 'dirichlet'
    
    #periodicity in x direction
    transf = Translation(w_domain*X)
    for domain, name in domain_dict.items():
        domain.faces.Min(X).name = name + '_x_min'
        domain.faces.Max(X).name = name + '_x_max'
        domain.faces.Min(X).Identify(domain.faces.Max(X), name+"_periodic_x", IdentificationType.PERIODIC, transf)

    #periodicity in y direction
    transf = Translation(w_domain*Y)
    for domain, name in domain_dict.items():
        domain.faces.Min(Y).name = name + '_y_min'
        domain.faces.Max(Y).name = name + '_y_max'
        domain.faces.Min(Y).Identify(domain.faces.Max(Y), name+"_periodic_y", IdentificationType.PERIODIC, transf)
    

    
    
    geo = OCCGeometry(Glue(domain_list))
    mesh = Mesh(geo.GenerateMesh(maxh=max(el_ag,el_air,el_sub)))

    return mesh