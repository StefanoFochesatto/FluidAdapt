from firedrake import *
from navierstokes import *
from animate import *
import os

os.chdir("/home/stefano/Desktop/AdaptiveFSI/Cylinder")

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
            DirichletBC(self.Z.sub(0), Constant((0.0, 0.0)), (14,)),  # circle
        ]

    def solve_step(self, replaceUold = False, updateT = False, write = False):
        """Solve one time step and write output."""
        solve(self.F == 0, self.up, bcs=self.bcs, nullspace=None, solver_parameters=self.sparams, options_prefix="")
        
        if updateT:
            self.t += self.dt
        if replaceUold:
            self.uold.interpolate(self.u)
        


    def get_hessian_metric(self):
        """Gets hessian metric of the current solution."""
        ux = self.u.sub(0)
        uy = self.u.sub(1)
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
    solver = AdaptiveNavierStokesSolver("cylinder.msh")
    solver.setup_problem()
    for i in range(100):

        # Initial Adaptation Step
        if i == 0:
            solver.solve_step(replaceUold = False, updateT = False)
            H = solver.get_hessian_metric()
            solver.mesh = adapt(solver.mesh, H)
            solver.setup_problem(prevMeshuold=solver.uold)
            solver.solve_step(replaceUold = True, updateT = False)
            solver.outfile.write(solver.u, solver.p, time=solver.t, adaptive = True)
            
        
        # Every buffer of 10 steps, solve the problem, avg hessian metric and update the mesh
        if i % 5 == 0:
            # Fill solution buffer
            storeU = solver.u
            metric_buffer = []
            for _ in range(5):
                solver.solve_step(replaceUold = True, updateT = False)
                metric = solver.get_hessian_metric()
                metric.normalise()
                metric_buffer.append(metric)
            # Compute Averaged metric
            H_avg = metric_buffer[0].copy(deepcopy=True)
            H_avg.average(*metric_buffer[1:],weights = [.4, .2, .2, .1, .1])
            H_avg.normalise()
            # Adapt mesh
            solver.mesh = adapt(solver.mesh, H_avg)
            # interpolate last solution to new adapted mesh
            solver.setup_problem(prevMeshuold=storeU)
            solver.outfile.write(solver.uold, solver.p, time=solver.t, adaptive = True)

            
            
        # resolve buffer on adapted mesh
        solver.solve_step(replaceUold = True, updateT = True)
        if i % 5 != 0:
            solver.outfile.write(solver.uold, solver.p, time=solver.t, adaptive = True)
        
        
        
    
    