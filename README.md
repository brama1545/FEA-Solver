# FEA Solver
 An FEA solving code that takes in an input file describing a static system and outputs deflections and internal forces based on finite element properties.
How to use this program:

This program models static FEA problems consisting of bar, spring, and/or beam elements.

Each element-type has its own class and stiffness matrices are defined for each.
An input file must be provided detailing the geometric coordinates and parameters of the system, as well as the element properties of each edge between two nodes.
Frames are combinations of beam and bar elements (Able to take axial, and bending loads), while grids are combos of beam and torsional spring elements (able to take bending moments and axial torques)
Loading and boundary conditions must be specified to fully define the system.

The program is designed to be as modular and flexible as possible. The only thing the user is responsible for is the input file and adding any distributed (non-nodal forces) to the workEq variable in the "solver.py" file using the work equivalence method. Additionally, output variables and changes to printing settings can be appended to the end of the "main" function in "solver.py" as desired.
Currently, only 1D elements can be analyzed. A second, independent program will soon be available which implements 2D planar elements such as constant strain triangles.

As the quarter progresses, new features will be incorporated (I hope I've allowed myself the flexibility to make changes quickly). Assuming everything else we learn this quarter follows the same process of building up a global stiffness matrix from a local stiffness matrix, and using {f} = [k]{d} to get the solution, this program should be a powerful way to solve a lot of different problems with different element types. Only the fileparser and edges files will need to be edited. The fileparsing needs to adapt as element types get more and more complicated and require more fields in the input file, and the edge function contains all the information as to how local stiffness matrices and transformation matrices are defined.
