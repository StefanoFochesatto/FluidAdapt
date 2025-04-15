from firedrake import *
from navierstokes import *
from animate import *
import os

os.chdir("/home/stefano/Desktop/FluidAdapt/MetricAdaptation/MetricAdvectionAdaptation")

class AdaptiveNavierStokesSolver:
    def __init__(self, mesh_file, dt=0.05, Re=400, H0=1.0, N=150, target_complexity=1500.0, h_min=1.0e-7, h_max=1.0):
        # Parameters for NS Problem
        self.mesh = Mesh(mesh_file)
        self.dt = dt
        self.Re = Re
        self.H0 = H0
        self.t = 0.0
        self.uold = None

        # Adaptation parameters
        self.target_complexity = target_complexity
        self.h_min = h_min
        self.h_max = h_max
        self.metric_buffer = []

        # Metric parameters
        self.metricparameters = {
            "dm_plex_metric": {
                "target_complexity": target_complexity,
                "p": 2.0,
                "h_min": h_min,
                "h_max": h_max,
                "boundary_tag": 14
            }
        }
        
        # Output file
        self.outfile = VTKFile("fluidTest.pvd")

    def metric_from_hessian(self, u):
        """Construct a hessian-based metric from a solution."""
        P1_ten = TensorFunctionSpace(self.mesh, "CG", 1)
        metric = RiemannianMetric(P1_ten)
        metric.set_parameters(self.metricparameters)
        metric.compute_hessian(u)
        metric.normalise()
        return metric

    def setup_problem(self, prevMeshuold = None):
        """Setup function spaces, forms, and boundary conditions."""
        self.Z, self.V, self.W, self.up = NSFunctions(self.mesh)
        self.u, self.p = self.up.subfunctions
        self.u.rename("u (velocity)")
        self.p.rename("p (pressure)")

        # Init uold,
        if self.uold is None:
            self.uold = Function(self.V).interpolate(as_vector([self.H0, 0.0]))
        elif prevMeshuold is not None:
            self.uold = Function(self.V).interpolate(prevMeshuold, allow_missing_dofs=True) # interpolate from previous mesh.
            
        
        self.F = NSTimeStepWeakForm(self.Z, self.up, self.uold, dt=self.dt, Re=self.Re)
        self.sparams = NSSolverParameters()

        # Boundary conditions
        ufar = Function(self.V).interpolate(as_vector([self.H0, 0.0]))
        self.bcs = [
            DirichletBC(self.Z.sub(0), ufar, (11, 13)),  # upstream, top, and bottom
            DirichletBC(self.Z.sub(0), Constant((0.0, 0.0)), (14)),  # circle
        ]

    def solve_step(self, replaceUold = False, updateT = False, write = False):
        """Solve one time step and write output."""
        solve(self.F == 0, self.up, bcs=self.bcs, nullspace=None, solver_parameters=self.sparams, options_prefix="")
        if write:
            self.outfile.write(self.u, self.p, time=self.t, adaptive = True)
        if updateT:
            self.t += self.dt
        if replaceUold:
            self.uold.interpolate(self.u)
        


    def get_hessian_metric(self, u = None):
        """Gets hessian metric of the current solution."""
        if u is None:
            ux = self.u.sub(0)
            uy = self.u.sub(1)
        else:
            ux = u.sub(0)
            uy = u.sub(1)
            
        Hx = self.metric_from_hessian(ux)
        Hy = self.metric_from_hessian(uy)
        Hx.normalise()
        Hy.normalise()
        H = Hx.copy(deepcopy=True)
        H.intersect(Hy)
        H.normalise
        return H
    
    
    
    def advect_metric(self, u_buffer):
        # function space for velocity
        V = VectorFunctionSpace(self.mesh, "CG", 2)
        # Constant function space for SUPG stabilized advection
        R = FunctionSpace(self.mesh, "R", 0)
        
        
        # Define boundary conditions for metric tensor
        P1_ten = TensorFunctionSpace(self.mesh, "CG", 1)
        h_max = self.h_max
        h_bc = Function(P1_ten).interpolate(Constant([[1.0 / h_max**2, 0.0], [0.0, 1.0 / h_max**2]]))
        
        ###### Setup for advection equation solve
        # Get initial metric tensor
        m0 = self.get_hessian_metric(u_buffer[0])
        m = RiemannianMetric(P1_ten) 
        metric_intersect = RiemannianMetric(P1_ten)
        metric_intersect.set_parameters(self.metricparameters)


        dt = Function(R).assign(self.dt)  # timestep size
        theta = Function(R).assign(0.5)  # Crank-Nicolson implicitness

        # SUPG stabilisation
        u0 = Function(V).interpolate(u_buffer[0])
        u = Function(V)
        
        D = Function(R).assign(0.1)
        h = CellSize(self.mesh)
        U = sqrt(dot(u, u))
        tau = 0.5 * h / U
        tau = min_value(tau, U * h / (6 * D))

        # Apply SUPG stabilisation
        phi = TestFunction(P1_ten)
        phi += tau * dot(u, grad(phi))

        # Variational form of the advection equation for the metric tensor
        trial = TrialFunction(P1_ten)
        a = inner(trial, phi) * dx + dt * theta * inner(dot(u, grad(trial)), phi) * dx
        L = inner(m0, phi) * dx - dt * (1 - theta) * inner(dot(u0, grad(m0)), phi) * dx
        bcs = DirichletBC(P1_ten, h_bc, (11, 12, 13))
        lvp = LinearVariationalProblem(a, L, m, bcs=bcs)
        lvs = LinearVariationalSolver(lvp)
        
        for ui in u_buffer[1:]:
            u.assign(ui)
            lvs.solve()
            m.enforce_spd(restrict_sizes=True, restrict_anisotropy=True)
            metric_intersect.intersect(m)
            m0.assign(m)
            u0.assign(u)
        
        metric_intersect.normalise()
        amesh = adapt(self.mesh, metric_intersect)
        return amesh
            
        
        
        

# The way this adaptation works is on each sub interval we will solve the problem, 
# use those velocity fields to advect the metric tensor on the same time interval, 
# intersect them all, and and adapt the mesh. 
# In order to do this i'll need a method which takes a buffer of velocity fields
# and advects the metric tensor on the same time interval.        
    

# Main execution
if __name__ == "__main__":
    solver = AdaptiveNavierStokesSolver("cylinder.msh")
    solver.setup_problem()
    for i in range(500):
        
        # Initial Adaptation Step
        if i == 0:
            solver.solve_step(replaceUold = False, updateT = False)
            H = solver.get_hessian_metric()
            solver.mesh = adapt(solver.mesh, H)
            solver.setup_problem(prevMeshuold=solver.uold)
            
    
        if i % 3 == 0:
            storeU = solver.uold
            u_buffer = []
            for _ in range(3):
                solver.solve_step(replaceUold = True, updateT = False)
                u_buffer.append(solver.u)
                
            # Advect Metric
            solver.mesh = solver.advect_metric(u_buffer)
            
            # interpolate last solution to new adapted mesh
            solver.setup_problem(prevMeshuold=storeU)
            
        # resolve buffer on adapted mesh
        solver.solve_step(replaceUold = True, updateT = True, write = True)
            

    
    