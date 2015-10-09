# dsmc
Implementation of a 3D Direct Simulation Monte Carlo code in Julia.

Status:

1. Import triangulated meshes from .ply ascii files.

2. Generate AMR grid around imported mesh (see two figures in doc)

3. Can push gas through simulation domain

4. Can output results to .vtk ascii files which can be read by Paraview.

To Do:
* Compute volume of the cut cells (cells that intersect with mesh of imported 3D model).
* All DSMC gas related capabilities, e.g. collisions between gas molecules, collisions with 3D surface, calculation of time step size etc.
* Speed up data dump to file.
* NAIF/SPICE integration
* Speed up code / make it parallel

