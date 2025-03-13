## specification

- DNS simulation of laminar channel flow
- Re = 200
- from uniform inflow to fully-developed laminar profile

## manipulation
3D mesh 50\*25\*1

Solve with `icoFoam`, the incompressible flow solver which gives transient results

increase  timestep `deltaT = 0.02`  to check stability condition

## operation

**OpenFOAM**

`blockMesh`

`simpleFoam`

**Paraview**

`paraFoam`, or open .foam file in Paraview GUI

plot *overline* to check velocity distribution

