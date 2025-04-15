from firedrake import *
from firedrake.output import VTKFile
import firedrake.cython.dmcommon as dmcommon
from firedrake.utils import IntType

from movement import *

import matplotlib.pyplot as plt
from firedrake.pyplot import triplot

import numpy as np


def filtermeshedges(mesh, indicator):
  
    # Create Section for DG0 indicator
    tdim = mesh.topological_dimension()
    entity_dofs = np.zeros(tdim+1, dtype=IntType)
    entity_dofs[:] = 0
    entity_dofs[-2] = 1 # -2 index for a 2d mesh is the edge stratum, 1 dof for selection matches Hdiv0T space 
    indicatorSect, _ = dmcommon.create_section(mesh, entity_dofs)

    # Pull Plex from mesh
    dm = mesh.topology_dm
    
    # Modify the rest of the code to use petsc mesh filtering
  
    # Create an adaptation label to mark cells for refinement
    dm.createLabel('filter')
    adaptLabel = dm.getLabel('filter')
    adaptLabel.setDefaultValue(0)

    # dmcommon provides a python binding for this operation of setting the label given an indicator function data array
    dmcommon.mark_points_with_function_array(
        dm, indicatorSect, 1, indicator.dat.data_with_halos, adaptLabel, 1)

    # Create a DMPlexTransform object to apply the refinement
    opts = PETSc.Options()

    opts['dm_plex_transform_active'] = 'filter'
    opts['dm_plex_transform_type'] = 'transform_filter' # 
    dmTransform = PETSc.DMPlexTransform().create(comm = mesh.comm)
    dmTransform.setDM(dm)
    # For now the only way to set the active label with petsc4py is with PETSc.Options() (DMPlexTransformSetActive() has no binding)
    dmTransform.setFromOptions()
    dmTransform.setUp()
    dmAdapt = dmTransform.apply(dm)
    
    # Labels are no longer needed, not sure if we need to call destroy on them. 
    dmAdapt.removeLabel('filter')
    dm.removeLabel('filter')
    dmTransform.destroy()
    
    # Remove labels to stop further distribution in mesh()
    # dm.distributeSetDefault(False) <- Matt's suggestion
    dmAdapt.removeLabel("pyop2_core")
    dmAdapt.removeLabel("pyop2_owned")
    dmAdapt.removeLabel("pyop2_ghost")
    # ^ Koki's suggestion

    # Pull distribution parameters from original dm
    distParams = mesh._distribution_parameters
    
    # Cast mesh as a 1D Mesh
    dmAdapt.setDimension(1)
    # Create a new mesh from the adapted dm
    refinedmesh = Mesh(dmAdapt, distribution_parameters = distParams, comm = mesh.comm)
    opts['dm_plex_transform_type'] = 'refine_regular'
    
    return refinedmesh






























def filtermesh(mesh, indicator):
  
    # Create Section for DG0 indicator
    tdim = mesh.topological_dimension()
    entity_dofs = np.zeros(tdim+1, dtype=IntType)
    entity_dofs[:] = 0
    entity_dofs[-1] = 1
    indicatorSect, _ = dmcommon.create_section(mesh, entity_dofs)

    # Pull Plex from mesh
    dm = mesh.topology_dm
    
    # Modify the rest of the code to use petsc mesh filtering
  
  
    # Create an adaptation label to mark cells for refinement
    dm.createLabel('filter')
    adaptLabel = dm.getLabel('filter')
    adaptLabel.setDefaultValue(0)

    # dmcommon provides a python binding for this operation of setting the label given an indicator function data array
    dmcommon.mark_points_with_function_array(
        dm, indicatorSect, 0, indicator.dat.data_with_halos, adaptLabel, 1)

    # Create a DMPlexTransform object to apply the refinement
    opts = PETSc.Options()

    opts['dm_plex_transform_active'] = 'filter'
    opts['dm_plex_transform_type'] = 'transform_filter' # 
    dmTransform = PETSc.DMPlexTransform().create(comm = mesh.comm)
    dmTransform.setDM(dm)
    # For now the only way to set the active label with petsc4py is with PETSc.Options() (DMPlexTransformSetActive() has no binding)
    dmTransform.setFromOptions()
    dmTransform.setUp()
    dmAdapt = dmTransform.apply(dm)
    
    # Labels are no longer needed, not sure if we need to call destroy on them. 
    dmAdapt.removeLabel('filter')
    dm.removeLabel('filter')
    dmTransform.destroy()
    
    # Remove labels to stop further distribution in mesh()
    # dm.distributeSetDefault(False) <- Matt's suggestion
    dmAdapt.removeLabel("pyop2_core")
    dmAdapt.removeLabel("pyop2_owned")
    dmAdapt.removeLabel("pyop2_ghost")
    # ^ Koki's suggestion

    # Pull distribution parameters from original dm
    distParams = mesh._distribution_parameters
    
    # Create a new mesh from the adapted dm
    refinedmesh = Mesh(dmAdapt, distribution_parameters = distParams, comm = mesh.comm)
    opts['dm_plex_transform_type'] = 'refine_regular'
    
    return refinedmesh


def plotBoundaryIDs(mesh):
    fig, axes = plt.subplots()
    firedrake.triplot(mesh, axes=axes)
    axes.legend()
    plt.show()



# Symmetric Gradient/Strain Tensor
def eps(u):
    return 0.5 * (grad(u) + grad(u).T)


# Isotropic hooke's law stress tensor
def sig(u, mu_, lambda_):
    return lambda_ * div(u) * Identity(2) + 2 * mu_ * eps(u)
    # equivalent:  lambda_ * tr(eps(u)) * Id + 2 * mu_ * eps(u)
    

# Cauchy Stress for computing stress tensor of fluid    
def sigF(u, p, mmu):
    return Constant(2) * mmu * eps(u) - p * Identity(2)     



def UpdateMesh(mesh, u):
    shape = mesh.coordinates.dat.data_with_halos.shape
    mesh.coordinates.dat.data_with_halos[:] += u.dat.data_with_halos.reshape(shape)
    mesh.bounding_box_cache = None
    UpdatePlex(mesh.topology_dm, mesh, u)
    

def UpdatePlex(plex,mesh, u):
    entity_dofs = np.zeros(mesh.topological_dimension() + 1, dtype=np.int32)
    entity_dofs[0] = mesh.geometric_dimension()
    coordinate_section = dmcommon.create_section(mesh, entity_dofs)[0]
    dm_coords = plex.getCoordinateDM()
    dm_coords.setDefaultSection(coordinate_section)
    local_coordinates_vec = dm_coords.createLocalVec()
    
    local_coordinates_vec.array[:] = np.reshape(
            mesh.coordinates.dat.data_with_halos,
            local_coordinates_vec.array.shape,
        )
    plex.setCoordinatesLocal(local_coordinates_vec)



if __name__ == "__main__":
    
    # Load the mesh with appropriate Surface IDs
    mesh = Mesh('fluidtunnelnew.msh')


    # Construct Indicator functions for fluid and solid
    V_W = FunctionSpace(mesh, "CG", 1)
    V_B = VectorFunctionSpace(mesh, "CG", 1)
    V_DG0_W = FunctionSpace(mesh, "DG", 0)
    V_DG0_B = FunctionSpace(mesh, "DG", 0)

    # Indicator function for the fluid
    I_W = Function(V_DG0_W)
    par_loop(("{[i] : 0 <= i < f.dofs}", "f[i, 0] = 1.0"),
            dx(1),
            {"f": (I_W, WRITE)})
    I_cg_W = Function(V_W)
    par_loop(("{[i] : 0 <= i < A.dofs}", "A[i, 0] = fmax(A[i, 0], B[0, 0])"),
            dx,
            {"A": (I_cg_W, RW), "B": (I_W, READ)})

    # Indicator functions for the solid
    I_B = Function(V_DG0_B)
    par_loop(("{[i] : 0 <= i < f.dofs}", "f[i, 0] = 1.0"),
            dx(2),

            {"f": (I_B, WRITE)})
    I_cg_B = Function(V_B)
    par_loop(("{[i, j] : 0 <= i < A.dofs and 0 <= j < 2}", "A[i, j] = fmax(A[i, j], B[0, 0])"),
            dx,
            {"A": (I_cg_B, RW), "B": (I_B, READ)})


    # Filter meshes 
    fluidmesh = filtermesh(mesh, I_W)
    solidmesh = filtermesh(mesh, I_B)


    # Mark boundaries for fluid
    # Create indicator functions for each boundary
    V = FunctionSpace(fluidmesh, "HDivT", 0)
    fluidindicators = []
    fluidIDs = [3, 4, 5, 6, 7, 9, 10, 11]
    x = SpatialCoordinate(fluidmesh)

    # Bottom boundaries (we'll handle them separately to avoid overlap)
    fluidindicators.append(Function(V).interpolate(conditional(abs(x[0]) < 1e-10, 1, 0))) # inlet
    fluidindicators.append(Function(V).interpolate(conditional(abs(x[0] - 4.0) < 1e-10, 1, 0))) # outlet
    fluidindicators.append(Function(V).interpolate(conditional(abs(x[1] - 1.0) < 1e-10, 1, 0))) # top
    fluidindicators.append(Function(V).interpolate(conditional(And(abs(x[1]) < 1e-10, x[0] <= 1.4), 1, 0))) # Bottom left
    fluidindicators.append(Function(V).interpolate(conditional(And(abs(x[1]) < 1e-10, x[0] >= 1.6), 1, 0)))



    # Fluid-solid boundary (combine solid right, top, left)
    solid_right = conditional(And(abs(x[0] - 1.6) < 1e-10, x[1] <= 0.5), 1, 0)
    solid_top = conditional(And(And(abs(x[1] - 0.5) < 1e-10, x[0] >= 1.4), x[0] <= 1.6), 1, 0)
    solid_left = conditional(And(abs(x[0] - 1.4) < 1e-10, x[1] <= 0.5), 1, 0)
    
    fluidindicators.append(Function(V).interpolate(solid_right))
    fluidindicators.append(Function(V).interpolate(solid_top))
    fluidindicators.append(Function(V).interpolate(solid_left))
    
    V = FunctionSpace(solidmesh, "HDivT", 0)
    solidindicators = []
    solidIDs = [8, 9, 10, 11]
    x = SpatialCoordinate(solidmesh)
    
    # Fluid-solid boundary (combine solid right, top, left)s
    solid_right = conditional(And(abs(x[0] - 1.6) < 1e-10, x[1] <= 0.5), 1, 0)
    solid_top = conditional(And(And(abs(x[1] - 0.5) < 1e-10, x[0] >= 1.4), x[0] <= 1.6), 1, 0)
    solid_left = conditional(And(abs(x[0] - 1.4) < 1e-10, x[1] <= 0.5), 1, 0)
    
    solidindicators.append(Function(V).interpolate(conditional(And(And(abs(x[1]) < 1e-10, x[0] >= 1.4), x[0] <= 1.6), 1, 0)))
    solidindicators.append(Function(V).interpolate(solid_right))
    solidindicators.append(Function(V).interpolate(solid_top))
    solidindicators.append(Function(V).interpolate(solid_left))

    b = Function(V).interpolate(solid_right + solid_top + solid_left)
    edgemesh = filtermeshedges(solidmesh, b)
    
    # BDids
    # Fluid Domain
    # 3 = inlet
    # 4 = outlet
    # 5 = upper
    # 6 = lowerleft
    # 7 = lowerright
    # Solid Domain
    # 8 = bottom
    # FSI
    # 11 = left
    # 10 = top
    # 9 = right
    fluidmesh = RelabeledMesh(fluidmesh, fluidindicators, fluidIDs)
    solidmesh = RelabeledMesh(solidmesh, solidindicators, solidIDs)
    plotBoundaryIDs(fluidmesh)
    plotBoundaryIDs(solidmesh)
    plt.show()
    
    
    
    outfileFluid = VTKFile("fluid.pvd", adaptive = True)
    outfileSolid = VTKFile("solid.pvd")
    
    t = 0
    Tmax = 8
    dt = 0.04
    uf0 = None
    ups0 = None
    um = None
    reffluidmesh = fluidmesh
    
    # Displacement space for mesh problem
    UM = VectorFunctionSpace(reffluidmesh, 'CG', 1)
    
    # Reference spaces for coupling
    RefUF = VectorFunctionSpace(reffluidmesh, 'CG', 1)
    RefPF = FunctionSpace(reffluidmesh, 'CG', 1)
    
    
    while t <= Tmax:
        # Setup for Fluid Problem
        if True:
            # Define Timestep size and Reynolds number
            rho = Constant(1.0)               # Fluid density [kg/m³]
            mu = Constant(.002)                  # Dynamic viscosity [Pa s]
            
            # Setup mixed function space for velocity and pressure
            UF = VectorFunctionSpace(fluidmesh, "CG", 2) # Velocity
            PF = FunctionSpace(fluidmesh, "CG", 1) # Pressure
            ZF = UF * PF
            upf = Function(ZF)
            
            if uf0 is None:
                uf0 = Function(UF)
            else:
                # Transfer previous velocity to new mesh
                uftransfer = Function(UF)
                uftransfer.dat.data_wo_with_halos[:] = uf0.dat.data_wo_with_halos[:]
                uf0 = uftransfer
                
            # Setup mesh Displacement
            if um is None:
                um = Function(UF)
                umv = Function(UF).interpolate(um/dt)
                
            else:
                # Transfer displacement computed on ref fluid mesh to new fluid mesh
                um_transfer = Function(VectorFunctionSpace(fluidmesh, 'CG', 1))
                um_transfer.dat.data_wo_with_halos[:] = um.dat.data_wo_with_halos[:]
                # Setup mesh Velocity
                umv = Function(UF).interpolate(um/dt)
            
            
            uf, pf = split(upf)
            vf, qf = TestFunctions(ZF)
            F = (
                rho * dot(uf, vf) * dx                     #  time derivative 
                + dt * mu * inner(grad(uf), grad(vf)) * dx                 # Viscous term (unchanged)
                + dt * rho * dot(dot(grad(uf), uf - umv), vf) * dx  # ALE convective term (ρ (u - w)·∇u)
                - dt * pf * div(vf) * dx                                   # Pressure term (unchanged)
                - rho * dot(uf0, vf) * dx                                  # Previous time step (unchanged)
                - div(uf) * qf * dx                                         # Continuity equation (unchanged)
            )
            
            bcs = [
                DirichletBC(ZF.sub(1), Constant(.5), (3,)),  # inlet
                #DirichletBC(ZF.sub(0), Constant((1.0, 0.0)), (3)), # sides
                DirichletBC(ZF.sub(0), Constant((0.0, 0.0)), (5, 6, 7)), # sides
                DirichletBC(ZF.sub(1), Constant(0.0), (4))  # outlet
            ]
            
            sparams = {
                    "snes_type": "vinewtonrsls",
                    "snes_rtol": 1.0e-6,
                    "snes_atol": 1.0e-9,
                    "ksp_type": "preonly",
                    "pc_type": "lu",
                    "pc_factor_mat_solver_type": "mumps",
                }
                
            uf, pf = upf.subfunctions
            
            uf.rename("uf (fluid velocity)")
            pf.rename("pf (fluid pressure)")
  
            print(f't = {t:.3f}:')
            solve(F == 0, upf, bcs=bcs, solver_parameters=sparams, options_prefix="")
            t += dt
            uf0.interpolate(uf)
            
            



        # at each time step we need to couple the velocity from the fluid to the solid through BCs

        # Setup for Solid Problem
        if True:
            # Setup material parameters and body forces. 
            rho = Constant(.25)*Constant(15)  # density
            g = Constant(0)  # gravity
            f = as_vector([0, -rho * g])  # body force
            mu_ = Constant(.25)*Constant(75)  
            lambda_ = Constant(.25)*Constant(125) 
            
            
            
            # Setup Mixed Function Space for displacement and velocity. 
            US = VectorFunctionSpace(solidmesh, "CG", 1) # Displacement
            PS = VectorFunctionSpace(solidmesh, "CG", 1) # Velocity
            ZS = US * PS
            ups = Function(ZS)
            
            if ups0 is None:
                ups0 = Function(ZS)
            else:
                ups0.interpolate(ups)
                
            us, ps = split(ups)
            us0, ps0 = split(ups0)
            vs, qs = TestFunctions(ZS)
            


            ####### Fully implicit coupled nonlinear weakform for displacement and velocity #######
            ### Main parts of weak form
            # rho(p') - div(sigma(u)) - f = 0
            F = rho * inner((ps - ps0), vs) * dx - dt * inner(sig(us, mu_, lambda_), eps(vs)) * dx - dt * inner(f, vs) * dx 
            # p - u' = 0
            F += inner((us - us0), qs) * dx  - dt * inner(ps, qs) * dx 
            
            
            ###### We need to define Neumann boundary conditions via the Piola transformed traction ###########
            # First thing we need to do is interpolate the fluid velocity into P1 space to match the solid problem. 
            UFP1 = VectorFunctionSpace(fluidmesh, "CG", 1)
            # Compatable fluid velocity. 
            ufP1 = Function(UFP1).interpolate(uf)

            # Compute the stress tensor on the transformed mesh
            sigF = Constant(2) * mu * eps(ufP1) - pf * Identity(2)
            TF = TensorFunctionSpace(fluidmesh, "CG", 1, (2, 2))
            sigFeval = Function(TF).interpolate(sigF)
            
            # Map Fluid Stress tensor back to reference mesh
            RefTF = TensorFunctionSpace(reffluidmesh, "CG", 1, (2, 2))
            sigFref = Function(RefTF)
            sigFref.dat.data_wo_with_halos[:] = sigFeval.dat.data_wo_with_halos[:]
            
            
            # Relevant terms for PiolaMap which transforms traction from transformed fluid to reference solid mesh        
            def DeformationGradient(u):
                return Identity(2) + grad(u)

            # We will want to pass stress tensor from the fluid, and the mesh velocities for the mesh problem, to compute J and F from umesh (29.15) 
            def PiolaMap(A, u):
                F = DeformationGradient(u)
                J = det(F)
                B = J*A*inv(F).T
                return B
                        
            # We compute the piola map on the reference mesh
            TFluid = PiolaMap(sigFref, um)
            
            # Define Neumann boundary condition            # Evaluate traction and transfer to solid mesh

            # Transfer solution from fluid to solid mesh along the interface
            TS = TensorFunctionSpace(solidmesh, "CG", 1, (2, 2))
            TE = TensorFunctionSpace(edgemesh, "CG", 1, (2, 2))
            # Interpolation exact since nodes are shared between meshes
            PiolaMapSigmaF = Function(RefTF).interpolate(TFluid) # on fluid mesh
            PiolaMapSigmaFedge = Function(TE).interpolate(PiolaMapSigmaF) # on FSI boundary (not really needed)
            TSolid= Function(TS).interpolate(PiolaMapSigmaFedge, allow_missing_dofs=True) # on solid mesh
            
            # Interpolate velocity from fluid to solid mesh
            ufP1ref = Function(RefUF)
            ufP1ref.dat.data_wo_with_halos[:] = ufP1.dat.data_wo_with_halos[:]
            VSolid = Function(US).interpolate(ufP1ref, allow_missing_dofs=True) 
            
            
            n = FacetNormal(solidmesh)
            F += inner(dot(TSolid, n), vs) * ds(9)  
            F += inner(dot(TSolid, n), vs) * ds(10)   
            F += inner(dot(TSolid, n), vs) * ds(11) 
            # Weakly enforced velocity dirchilet bc
            #F += inner(ps - VSolid, qs) * ds(9)
            #F += inner(ps - VSolid, qs) * ds(10)
            #F += inner(ps - VSolid, qs) * ds(11)
            
            
            
            bcs = [DirichletBC(ZS.sub(0), Constant([0, 0]), 8),  # Fix solid base
                DirichletBC(ZS.sub(1), VSolid, 9),  # Fix solid velocities on FSI
                DirichletBC(ZS.sub(1), VSolid, 10), 
                DirichletBC(ZS.sub(1), VSolid, 11)  
                ] 
            
            sp = {
                "snes_type": "newtonls",
                "ksp_type": "preonly",
                "pc_type": "lu",
                "snes_converged_reason": None,
                "snes_monitor": None,
                "ksp_monitor": None,
            }
            us, ps = ups.subfunctions
            solve(F == 0, ups, bcs=bcs, solver_parameters=sp)
            
            us.rename("u solid displacement")
            ps.rename("p solid velocity")
            TSolid.rename("Traction Boundary Condition")
            VSolid.rename("Velocity Boundary Condition")
            outfileSolid.write(us, ps,VSolid, TSolid, t = t - dt)
        
                    
            
            
        # The fluid mesh needs to be updated with the displacement from the solid. The mesh problem needs to be formulated on the reference fluid mesh
        if True: 
            mu_M = Constant(3.8461)
            lambda_M = Constant(5.76)
            us, ps = ups.subfunctions
            umBDRY = Function(UM, name = "Mesh Dirchilet Boundary Conditions").interpolate(us, allow_missing_dofs=True) #interpolate second order displacement to first order mesh space 
            
            def SigmaM(um):
                return mu_M * (grad(um) + grad(um).T) + lambda_M * tr(grad(um)) * Identity(2)
            
            um0 = Function(UM, name = "Mesh Displacement")
            um = Function(UM, name = "Mesh Displacement")
            vm = TestFunction(UM)
            
            F = inner(um, vm)*dx  - .5*dt*inner(SigmaM(um), grad(vm) + grad(vm).T)*dx
            F += -inner(um0, vm)*dx + .5*dt*inner(SigmaM(um0), grad(vm) + grad(vm).T)*dx
            
            bcs = [
                DirichletBC(UM, Constant((0.0, 0.0)), [3, 4, 5, 6, 7]), 
                DirichletBC(UM, umBDRY, [9, 10, 11])
            ]
                        
            sp = {
                "snes_type": "newtonls",
                "ksp_type": "preonly",
                "pc_type": "lu",
                "snes_converged_reason": None,
                "snes_monitor": None,
                "ksp_monitor": None,
            }
            solve(F == 0, um, bcs=bcs, solver_parameters=sp)
            um0.assign(um)
            fluidmesh = Mesh(reffluidmesh.coordinates.copy(deepcopy=True))
            UpdateMesh(fluidmesh, um)
            
            test = Function(VectorFunctionSpace(fluidmesh, "CG", 2))
            test.dat.data[:] = uf.dat.data[:]
            test.rename("velocitu")
            outfileFluid.write(test, t = t, adaptive = True)
            
            
            print("test")
            
            
            
            
            
            # umesh = Function(RefV).interpolate(us, allow_missing_dofs=True) #dont' know why I have to multiply by 100
            # umesh.rename("mesh displacement")
            # VTKFile("meshprobBDRYCond.pvd").write(umesh)

            # mover = SpringMover(reffluidmesh, .05 ,method = 'lineal')
            # fsi = Function(mover.coord_space)
            # moving_boundary = DirichletBC(mover.coord_space, fsi, [9, 10, 11])


            # def update_boundary_displacement(t):
            #     coord_data = mover.mesh.coordinates.dat.data
            #     for i in moving_boundary.nodes:
            #         x, y = coord_data[i]
            #         fsi.dat.data[i] = umesh.at(x, y)
            
            # fixed_boundaries = DirichletBC(mover.coord_space, 0, [3, 4, 5, 6, 7])
            # boundary_conditions = (fixed_boundaries, moving_boundary)
            
            # mover.move(t - dt, update_boundary_displacement=update_boundary_displacement, boundary_conditions=boundary_conditions,)
            
            # um = Function(RefV)
            # um.dat.data_wo_with_halos[:] = mover.mesh.coordinates.dat.data_wo_with_halos[:] - fluidmesh.coordinates.dat.data_wo_with_halos[:]
            
            # fluidmesh = mover.mesh
            # VTKFile('newmesh.pvd').write(fluidmesh)
            
            