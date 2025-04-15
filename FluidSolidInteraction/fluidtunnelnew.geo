// fluidtunnel.geo

// Define characteristic length (global, for simplicity)
lc = 0.3;

// Define points
Point(1) = {0, 0, 0, lc};   // Bottom left corner of the domain
Point(2) = {4, 0, 0, lc};   // Bottom right corner of the domain
Point(3) = {4, 1, 0, lc};   // Top right corner of the domain
Point(4) = {0, 1, 0, lc};   // Top left corner of the domain
Point(5) = {1.4, 0, 0, lc}; // Bottom left corner of the solid
Point(6) = {1.6, 0, 0, lc}; // Bottom right corner of the solid
Point(7) = {1.6, 0.5, 0, lc}; // Top right corner of the solid
Point(8) = {1.4, 0.5, 0, lc}; // Top left corner of the solid


// Define lines
Line(1) = {1, 5};  // Bottom of the domain to the left of the solid
Line(2) = {6, 2};  // Bottom of the domain to the right of the solid

Line(10) = {1, 2}; // Bottom of the entire domain
Line(3) = {2, 3};   // Right side of the domain
Line(4) = {3, 4};   // Top of the entire domain
Line(5) = {4, 1};   // Left side of the domain

Line(6) = {5, 6};   // Bottom of the solid - Important.
Line(7) = {6, 7};   // Right side of the solid
Line(8) = {7, 8};   // Top of the solid
Line(9) = {8, 5};   // Left side of the solid


// Define curve loops
Curve Loop(1) = {1, 6, 2, 3, 4, 5};   // Outer boundary (fluid domain)
Curve Loop(2) = {6, 7, 8, 9};   // Inner boundary (solid)


// Define surfaces
Plane Surface(1) = {1,-2}; // Fluid, the `-2` subtracts the solid region's boundary (Curve Loop 2)
Plane Surface(2) = {2};       // Solid,  this refers to the *inside* of  Curve Loop 2.  This surface definition is needed for solid

// Define Physical Surfaces
Physical Surface("fluid") = {1};   // Fluid (includes the entire domain, with the solid's area excluded)
Physical Surface("solid") = {2};   // Solid (the area inside the solid's boundary)


// Physical lines (for boundary conditions)
// Physical Line("bottom") = {10}; // Bottom of the domain
// Physical Line("outlet") = {3};  // Right side of the domain
// Physical Line("inlet") = {5};   // Left side of the domain
// Physical Line("top") = {4};     // Top of the domain

// Physical Line("bottom_left") = {1};  // Bottom of the domain to the left of the solid
// Physical Line("bottom_right") = {2}; // Bottom of the domain to the right of the solid
// Physical Line("solid_bottom") = {6}; // Bottom of the solid

// Physical Line("fluid_solid_boundary") = {7, 8, 9}; // the sides and top of the solid, excluding the bottom which is a part of 'bottom'

Mesh.CharacteristicLengthFactor = 0.3;