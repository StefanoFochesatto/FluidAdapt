## Notes

This will be an implementation of Chapter 29 adaptive FSI solver from the FeniCS book.

In an FSI problem, physical quantities between fluid and solids are transferred along a FS boundary traction forces are given by normal stresses at the FS boundary, 
$$\sigma_F \cdot n_F = -\sigma_S \cdot n_S $$
This is essentially saying that the force on the fluid is equal and opposite to the force on the solid.

We also require that the velocity of the fluid and solid are equal at the FS boundary,
$$u_F = u_S$$
This is called kinematic continuity. 

### Lagrangian Framework for Structural Mechanics
In structural mechanics we care about the displacement field, in most simulations this is 
represented in a Lagrangian framework, i.e we are tracking nodes as they are moved. 
Let $x_0$ be a material point
let $\phi(x_0, t)$ be the position of the material point at time $t$
Then displacement the displacement of $x_0$ is given by $u_s(x_0, t) = \phi(x_0, t) - x_0$
Naturally, can differentiate and we get a non-singular Jacobi matrix $f = \nabla \phi$ and jacobi determinant $j = det f$

#### Compressible St. Venant-Kirchhoff Model
Chapter 29 focuses on the compressible St. Venant-Kirchhoff model where stress is described by the energy functional, 
$$ \psi(e) = \mu_s tr(e^2) + \frac{\lambda_s}{2} (tr(e))^2 $$
where $e = \frac{1}{2}(f^T f - I)$ is the Green-Lagrange strain tensor.  The Piola-Kirchhoff stress tensor is then given by
$$ \sigma_s(u_s) = f \cdot \frac{\partial \psi}{\partial e} =f\cdot( 2 \mu_s e + \lambda_s tr(e) I)$$

Conservation laws for the compressible St. Venant-Kirchhoff model are given by
$$ \rho_s \ddot{u_s} - \text{div} (\sigma_s(u_s)) = b_s$$
$$\dot{\rho_s} = 0$$

The text states that mass conservation is usually omitted in the compressible case in the Lagrangian formulation.

### Eulerian Framework for Fluid Mechanics

In fluid mechanics we care about fluid velocity $u_f$ and pressure $p_F$. These fields are naturally posed in the Eulerian framework where he motion of a body is related to a fixed spatial point $x$ as 
$$u_f(x, t) = u_f(\phi(x_0), t)$$ 
Thus the material time derivative of a field $y$ in the Eulerian framework is given by
$$d_t(y) = \dot{y} + \underbrace{\text{grad} (y) \cdot u_f}_{\text{convective term}}$$
The convective term is the term that makes the Eulerian framework difficult to work with from a material point of view, as opposed to the Lagrangian framework

#### Incompressable Newtonian Fluids Model

In the Newtonian fluid model, the Cauchy stress tensor $\sigma_f$ is given by
$$\sigma_f(u_f, p_f) = 2\mu_f\epsilon(u_f) - p_fI$$

Where $\mu_f$ is the dynamic viscosity and $\epsilon(\cdot)$ is the symmetric gradient. The fluid is then described by the incompressable Navier-Stokes equations
$$\rho_f((\dot u_f) + \text{grad}(u_f) \cdot u_f) - \text{div} (\sigma_f(u_f, p_f)) = b_f$$
$$\text{div}(u_f) = 0$$

### FSI and the ALE Framework

To bring both frameworks together there are two key steps. At the FS interface the force from the fluid is transferred to the the solid via the Piola Map 
$$ (j \sigma_f \cdot f^{-T})\cdot n_f = -\sigma_s \cdot n_s$$

The structure in the material domain must be tracked in the spatial fluid domain, and therefore the spatial fluid mesh needs to update, in an intelligent way which does not destroy the quality of the mesh. The second key step is the mesh smoothing procedure which solves the aformentioned problem. 

##### An Arbitrary Reference Frame
Let $\Omega$ be the reference(undeformed) computational domain and let $\Omega$ be partitioned into two sets $\Omega_s$ and $\Omega_f$. Let $\omega(t)$ denote the current deforemed computational domain which is also partitioned into $\omega_s(t)$ and $\omega_f(t)$ for all $t \in [0, T]$. We'll define the FS boundary as $\Gamma_{fs}$ and $\gamma_{fs}(t)$ on the reference and current domains respectively. We'll use upper and lowercase letters to define respective fields on the reference and current domains. We'll consider the map for mapping a point $X$ in the reference domain to a point $x$ in the current domain as
$$ \Phi: \Omega \rightarrow \omega(t)$$ 
which is given by 
$$ \Phi(X, t) \begin{cases}
\Phi_s(X, t) , \quad \forall X \in \Omega_s \\
\Phi_m(X, t) , \quad \forall X \in \Omega_f
\end{cases}
$$
$$ \Phi_s(X, t) = X + U_s(X, t) $$
$$ \Phi_m(X, t) = X + U_m(X, t) $$
Here $(U_s, U_m)$ are the solutions to the structure problem and the mesh problem which resolves the issues with needing to deform $\Omega_f$ over time as the solid deforms, to maintain a good mesh. 
There are many ways to formulate the mesh problem, the Movement library offers many options I've seen several presentations of ALE that use a concept called Laplacian Smoothing. Chapter 29 of the fenics book proposes we treat the fluid domain as a  linearly elastic solid (this is the lineal springs approach which is also implemented in Movement) with the stress tensor is given by, 
$$\Sigma_m(U_m) = \mu_m (\text{grad}(U_m) + \text{grad}(U_m)^T) + \gamma_m \text{tr}(\text{grad}(U_m))I$$









## GamePlan:

1. Write a fluid and solid solvers that can be used for each problem independently. 
2. Use the FSI tutorial in the firedrake docs to setup the mesh and BC conditions and FS boundary
3. Use the Movement Library to solve the mesh problem

Potential issues. We need to figure out what the FenicS book is saying when page 556 with respect to the convective term of the fluid problem when you solve it on the deformed domain. 
If we use the inverse map of $\Phi$ on  we can solve the entire problem on the reference domain

So here's how we're going to do this, we need to find a way to have gmsh spit out to meshes for each fluid and solid subdomain, and the entire reference domain, with all labels for the boundaries and FSI boundary. 

Then we will solve the ALE formulation of the Stokes Equation on the fluid domain, likely implicit solver with mesh velocity correction. 

Then we will use the FSI boundary condition to setup the solids problem,

Then we will use the movement library to solve the mesh problem.

In order to generate the separate meshes we will need to use the dmplextransformation functionality to filter the mesh for each subdomain. (this will be better for adaptivity we we decide decoupled adaptivity)
   


   The finite element approach to modeling time-dependent hyperelasticity in structural mechanics involves several key steps: spatial discretization, time discretization, linearization, and solving the resulting nonlinear system. Below, I outline the mathematical formulation and provide a **Firedrake** implementation for implicit time stepping.

---

## **1. Governing Equations for Hyperelasticity**
Hyperelastic materials are characterized by a strain energy density function \( \Psi(\mathbf{F}) \), where \( \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} \) is the deformation gradient, and \( \mathbf{u} \) is the displacement.

### **Strong Form (Momentum Balance)**

$$
\rho \ddot{\mathbf{u}} - \nabla \cdot \mathbf{P}(\mathbf{F}) = \mathbf{f}
$$
where:
- \rho  = density
- \( \mathbf{P} = \frac{\partial \Psi}{\partial \mathbf{F}} \) = First Piola-Kirchhoff stress
- \( \mathbf{f} \) = body force

### **Weak Form (Virtual Work Principle)**
Multiply by a test function \( \mathbf{v} \in V \) and integrate:
$$
\int_{\Omega_0} \rho \ddot{\mathbf{u}} \cdot \mathbf{v} \, dV + \int_{\Omega_0} \mathbf{P}(\mathbf{F}) : \nabla \mathbf{v} \, dV = \int_{\Omega_0} \mathbf{f} \cdot \mathbf{v} \, dV + \int_{\partial \Omega_0} \mathbf{t} \cdot \mathbf{v} \, dS
$$
where \( \mathbf{t} \) is the traction.

---

## **2. Time Discretization (Implicit Newmark-β)**
For implicit time stepping, we use the **Newmark-β** method:
$$
\ddot{\mathbf{u}}_{n+1} = \frac{\mathbf{u}_{n+1} - \mathbf{u}_n - \Delta t \dot{\mathbf{u}}_n}{\beta \Delta t^2} - \frac{1-2\beta}{2\beta} \ddot{\mathbf{u}}_n
$$
$$
\dot{\mathbf{u}}_{n+1} = \dot{\mathbf{u}}_n + \Delta t \left( (1-\gamma) \ddot{\mathbf{u}}_n + \gamma \ddot{\mathbf{u}}_{n+1} \right)
$$
where \( \beta = 0.25 \) (trapezoidal rule) and \( \gamma = 0.5 \) for unconditional stability.

Substitute \( \ddot{\mathbf{u}}_{n+1} \) into the weak form to get a nonlinear problem in \( \mathbf{u}_{n+1} \).

---

## **3. Firedrake Implementation**
Below is a **Firedrake** code snippet solving a time-dependent hyperelastic problem with implicit Newmark-β.

### **Key Steps:**
1. Define the mesh and function spaces.
2. Specify the hyperelastic material model (e.g., Neo-Hookean).
3. Apply time discretization (Newmark-β).
4. Solve the nonlinear problem at each time step.

### **Code Example:**
```python
from firedrake import *
import numpy as np

# Define mesh and function space
mesh = UnitCubeMesh(10, 10, 10)
V = VectorFunctionSpace(mesh, "CG", 1)

# Material parameters
rho = Constant(1.0)       # Density
mu = Constant(1.0)        # Shear modulus
lmbda = Constant(1.0)     # Lamé parameter

# Time stepping parameters
dt = 0.1
T = 1.0
beta = 0.25
gamma = 0.5

# Newmark-beta coefficients
a0 = Constant(1.0 / (beta * dt**2))
a1 = Constant(gamma / (beta * dt))
a2 = Constant(1.0 / (beta * dt))
a3 = Constant(1.0 / (2 * beta) - 1.0)

# Define trial and test functions
u = Function(V)       # Displacement at n+1
v = TestFunction(V)   # Test function

# Previous time step values
u_n = Function(V)     # Displacement at n
u_dot_n = Function(V) # Velocity at n
u_ddot_n = Function(V) # Acceleration at n

# Define the deformation gradient F = I + grad u
I = Identity(3)
F = I + grad(u)

# Neo-Hookean strain energy
J = det(F)
C = F.T * F
Ic = tr(C)
Psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# First Piola-Kirchhoff stress
P = diff(Psi, F)

# Inertia term (Newmark-beta discretized)
inertia = rho * a0 * inner(u - u_n, v) * dx
damping = rho * a1 * inner(u_dot_n, v) * dx
acceleration = rho * inner(a3 * u_ddot_n, v) * dx

# Internal and external forces
internal_force = inner(P, grad(v)) * dx
external_force = inner(Constant((0, 0, -0.1)), v) * dx  # Body force

# Weak form residual
residual = inertia + damping + acceleration + internal_force - external_force

# Boundary conditions (clamped at bottom)
bc = DirichletBC(V, Constant((0, 0, 0)), "bottom")

# Nonlinear solver setup
problem = NonlinearVariationalProblem(residual, u, bcs=bc)
solver = NonlinearVariationalSolver(problem)

# Time stepping loop
t = 0.0
while t < T:
    solver.solve()
    
    # Update acceleration, velocity, displacement
    u_ddot = a0*(u - u_n) - a2*u_dot_n - a3*u_ddot_n
    u_dot = u_dot_n + dt*((1-gamma)*u_ddot_n + gamma*u_ddot)
    
    # Update previous step
    u_n.assign(u)
    u_dot_n.assign(u_dot)
    u_ddot_n.assign(u_ddot)
    
    t += dt
    print(f"Time = {t:.3f}, Max displacement = {np.max(u.vector().array()):.3e}")
```

---

## **4. Key Takeaways**
- **Weak form** combines inertia, damping, and hyperelastic stress.
- **Implicit Newmark-β** ensures stability for large time steps.
- **Firedrake** automates finite element assembly and nonlinear solving.
- The **Neo-Hookean model** is used here, but other hyperelastic laws can be substituted.

### **Extension to Fluid-Structure Interaction (FSI)**
For FSI, you would:
1. Couple this structural solver with a Navier-Stokes solver (e.g., using **PETSc’s DMComposite**).
2. Apply interface conditions (traction continuity, no-slip velocity).

Would you like an example of a coupled FSI problem in Firedrake as well?