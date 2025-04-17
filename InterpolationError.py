from firedrake import *
from animate import *

mesh = RectangleMesh(nx=10, ny=10, Lx=1, Ly=1, originX = -1, originY = -1)

S = FunctionSpace(mesh, "CG", 1)
x, y = SpatialCoordinate(mesh)
u = Function(S).interpolate(sqrt(1 - x**2))

VTKFile("original.pvd").write(u)


h_min=1.0e-12
h_max=1



target_complexity = len(u.dat.data)-115
metricparameters = {
            "dm_plex_metric": {
                "target_complexity": target_complexity,
                "p": 2.0, 
                "h_min": h_min,
                "h_max": h_max
                }
        }


P1_ten = TensorFunctionSpace(mesh, "CG", 1)
metric = RiemannianMetric(P1_ten)
metric.set_parameters(metricparameters)
metric.compute_hessian(u)
metric.normalise()
meshAdapt = adapt(mesh, metric)

SAdapt = FunctionSpace(meshAdapt, "CG", 1)
x, y = SpatialCoordinate(meshAdapt)
uAdapt = Function(SAdapt).interpolate(sqrt(1 - x**2))
VTKFile("adapted.pvd").write(uAdapt)


