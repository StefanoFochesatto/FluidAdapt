from firedrake import *
from navierstokes import *
from animate import *
from movement import *
import os

os.chdir("/home/stefano/Desktop/FluidAdapt/MovementAdaptation/ClassicalAdaptation")

class AdaptiveNavierStokesSolver:
    def __init__(self, mesh_file, dt=0.05, Re=400, H0=1.0, N=150, target_complexity=100.0, h_min=1.0e-7, h_max=1.0):
        # Parameters for NS Problem
        self.buffermesh = mesh
        self.bufferu = None
        self.mesh = mesh
        self.dt = dt
        self.rtol = 1.0e-02
        self.Re = Re
        self.H0 = H0
        self.t = 0.0
        self.uold = None
        self.mu = None
        
        
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
                "h_max": h_max
            }
        }
        

        
        # Output file
        self.outfile = VTKFile("fluidTest.pvd")



    def setup_problem(self, prevMeshuold = None, mu = None):
        """Setup function spaces, forms, and boundary conditions."""
        self.Z, self.V, self.W, self.up = NSFunctions(self.mesh)
        self.u, self.p = self.up.subfunctions
        self.u.rename("u (velocity)")
        self.p.rename("p (pressure)")
        

        # Init uold,
        if self.uold is None:
            self.uold = Function(self.V).interpolate(as_vector([0.0, 0.0]))
        elif prevMeshuold is not None:
            self.uold = Function(self.V)
            self.uold.dat.data[:] = prevMeshuold.dat.data[:] #Transfer velocities
        
        # Init mesh velocity
        if mu is None:
            self.mu = Function(self.V).interpolate(as_vector([0.0, 0.0]))
        else:
            self.mu = Function(self.V).interpolate(mu) #mesh velocity will be a CG1 space, we interpolate to fit our form
        
        self.F = NSTimeStepWeakFormALE(self.Z, self.up, self.uold, self.mu ,dt=self.dt, Re=self.Re)
        self.sparams = NSSolverParameters()

        # Boundary conditions
        ufar = Function(self.V).interpolate(as_vector([self.H0, 0.0]))
        self.bcs = [
        DirichletBC(self.Z.sub(0), Constant((1.0, 0.0)), (4,)),  # top
        DirichletBC(self.Z.sub(0), Constant((0.0, 0.0)), (1, 2))  # sides
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
            
    
    def metric_from_hessian(self, u):
        """Construct a hessian-based metric from a solution."""
        P1_ten = TensorFunctionSpace(u._function_space.mesh(), "CG", 1)
        metric = RiemannianMetric(P1_ten)
        metric.set_parameters(self.metricparameters)
        metric.compute_hessian(u)
        metric.normalise()
        return metric

    
    
    
    
    
    
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
    
        
    
        
# Main execution
if __name__ == "__main__":
    mesh = UnitSquareMesh(10, 10)
    solver = AdaptiveNavierStokesSolver(mesh)
    solver.setup_problem()
    
    
        
    P1_ten = TensorFunctionSpace(solver.mesh, "CG", 1)
    H = RiemannianMetric(P1_ten)
    Q = FunctionSpace(solver.mesh, "CG", 1)
    Hnorm = Function(Q, name="Hnorm")
    alpha = Constant(10.0)
    
    
    
    
    def monitor(mesh):
        # Project the stashed solution data onto the current mesh
        Utmp = VectorFunctionSpace(solver.buffermesh, "CG", 2)
        solver.bufferu = Function(Utmp)
        solver.bufferu.dat.data[:] = solver.u.dat.data
        solver.u.project(solver.bufferu)
            

        # Compute the monitor expression using the Hessian
        H = solver.get_hessian_metric(solver.u)
        Hnorm.interpolate(sqrt(inner(H, H)))
        Hnorm_max = Hnorm.dat.data.max()
        m = 1 + alpha * Hnorm / Hnorm_max
        
        solver.buffermesh = Mesh(mesh.coordinates.copy(deepcopy=True))
        VTKFile('test.pvd').write(mesh)
        return m
    
    mover = MongeAmpereMover(solver.mesh, monitor, method="quasi_newton", rtol=solver.rtol)
    for i in range(500):
        if i == 0:
            solver.solve_step(replaceUold = False, updateT = False, write=False)
            solver.bufferu = solver.u
            solver.buffermesh = solver.mesh
            storeu = solver.u
            # update mesh
            mover.move()
            solver.mesh = mover.mesh
            
        solver.solve_step(replaceUold = True, updateT = True, write=True)
        solver.bufferu = solver.u
        solver.buffermesh = solver.mesh
        storeu = solver.u
        # update mesh
        mover.move()
        # Compute displacment
        MU = VectorFunctionSpace(solver.mesh, "CG", 1)
        disp = Function(MU)
        
        disp.dat.data[:] = mover.mesh.coordinates.dat.data[:] - solver.mesh.coordinates.dat.data[:]
        mu = disp/solver.dt
        
        solver.mesh = mover.mesh
        solver.setup_problem(storeu, mu)
        

        
    
    