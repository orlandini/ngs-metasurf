def drawDomainsAll(mesh, drawFunc=None, block=True):
    drawDomains(mesh, name=set(mesh.GetMaterials()), drawFunc=drawFunc, block=block) 
def drawDomains(mesh, name=None, drawFunc=None, block=True):
    from ngsolve import CF, Integrate

    if drawFunc == None:
        from ngsolve import Draw
        drawFunc = Draw

    if name==None:
        vals = dict(enumerate(mesh.GetMaterials()))
        vals = {v: k for k, v in vals.items()}
        drawFunc(mesh.MaterialCF(vals), mesh, "domains")
    else:
        if not hasattr(name, "__iter__") or type(name) == str:
            name = [name]
        
        for n in name:

            domCF = mesh.MaterialCF({n:1}, default=0)
            print(f"---- {n} ----- {Integrate(domCF,mesh)} ")
            drawFunc(domCF, mesh, "domains")
            
            if block:
                input("press Enter")
            
def drawBndAll(mesh, block=True, old=None, drawFunc=None, useBBnd=False, **kwargs):
    [drawBnd(mesh, x, block, old, drawFunc=drawFunc, useBBnd=useBBnd, **kwargs) for x in set(mesh.GetBoundaries() if not useBBnd else mesh.GetBBoundaries())]

def drawBnd(mesh, name="bottom|right|top|left|ibot|itop|interface|ileft|iright|natural", block=False, old=None, drawFunc=None, useBBnd=False, **kwargs):
    from ngsolve import CF, H1, GridFunction, BND, BBND, Integrate

    if old == None:
        old = True if mesh.dim <= 2 else False
    if drawFunc == None:
        from ngsolve import Draw
        drawFunc = Draw

    val = {bnd:1 for bnd in name.split("|")}

    if old:
        if useBBnd:
            fes = H1(mesh, dirichlet_bbnd=name)
            sol = GridFunction(fes, "bbnd")
            sol.Set(mesh.BBoundaryCF({name:1}, default=0), VOL_or_BND=BBND)
        else:
            fes = H1(mesh, dirichlet=name)
            sol = GridFunction(fes, "bnd")
            sol.Set(mesh.BoundaryCF({name:1}, default=0), VOL_or_BND=BND)
        drawFunc(sol, **kwargs)
        print("-----", name, sum(sol.vec))

    else:
        if useBBnd:
            bnd_CF = mesh.BBoundaryCF({name:1}, default=0)
        else:
            bnd_CF = mesh.BoundaryCF({name:1}, default=0)
        drawFunc(bnd_CF, mesh, "bnd", draw_vol=False)
        print("-----", name, Integrate(bnd_CF, mesh, VOL_or_BND=(BND if not useBBnd else BBND)))

    if block:

        input("press Enter")
