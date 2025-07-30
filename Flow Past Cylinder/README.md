flow past a cylinder


blockMesh
checkMesh
decomposePar -force
mpirun -n 4 icoFoam -parallel
reconstructPar -latestTime
