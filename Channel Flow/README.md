## specification

- DNS (laminar) simulation of laminar channel flow
- case preset Re = 200
- from uniform inflow to fully-developed laminar profile

## manipulation
3D mesh 50\*25\*1

Solve with `icoFoam`, the incompressible (laminar) flow solver for transient results

Increase  timestep `deltaT = 0.02`  to check stability condition

## command

**OpenFOAM**

`blockMesh`

`simpleFoam`

**Paraview**

`paraFoam`, or open .foam file in Paraview GUI

plot *overline* to check velocity distribution

