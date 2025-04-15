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

1. Write a fluid and solid and mesh problem solvers that can be used for each problem independently. 
2. Use the FSI tutorial in the firedrake docs to setup the mesh and BC conditions and FS boundary