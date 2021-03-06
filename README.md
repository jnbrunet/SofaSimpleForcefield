# Simple implementation of a Saint-Venant-Kirchhoff force field

This plugin can be used as a simple demonstration of how SOFA 
"force fields" components are built up. In this case, the force
field API will provide the internal potential energy,
residual and jacobian of a Saint-Venant-Kirchhoff material assembled
on linear tetrahedral meshes. These three function will be used 
by the system assembler and ODE solvers of SOFA.

## Compiling with Ubuntu LTS 20.04
### SOFA (skip this if you already compiled SOFA)

```console
user@host:~$ sudo apt install qtbase5-dev libqt5charts5-dev libqt5opengl5-dev libopengl0 libeigen3-dev libglew-dev zlib1g-dev libboost-dev libboost-filesystem-dev g++ cmake git
user@host:~$ export SOFA_SRC=/opt/sofa_src
user@host:~$ export SOFA_ROOT=/opt/sofa
user@host:~$ git clone https://github.com/sofa-framework/sofa.git $SOFA_SRC
user@host:~$ cmake -DCMAKE_INSTALL_PREFIX=$SOFA_ROOT -DCMAKE_BUILD_TYPE=Release -S $SOFA_SRC -B $SOFA_SRC/build
user@host:~$ cmake --build $SOFA_SRC/build -j4
user@host:~$ cmake --install $SOFA_SRC/build
```

### SofaSimpleForcefield
```console
user@host:~$ export SSFF_SRC=/opt/SofaSimpleForceField_src
user@host:~$ git clone https://github.com/jnbrunet/SofaSimpleForcefield.git $SSFF_SRC
user@host:~$ cmake -DCMAKE_PREFIX_PATH=$SOFA_ROOT/lib/cmake -DCMAKE_INSTALL_PREFIX=$SOFA_ROOT/plugins/SofaSimpleForceField -DCMAKE_BUILD_TYPE=Release -S $SSFF_SRC -B $SSFF_SRC/build
user@host:~$ cmake --build $SSFF_SRC/build -j4
user@host:~$ cmake --install $SSFF_SRC/build
```

## Running the cantilever example scene
```console
user@host:~$ $SOFA_ROOT/bin/runSofa -l SofaSimpleForceField $SSFF_SRC/cantilever_beam.scn
```

## Result
![image](https://user-images.githubusercontent.com/6951981/127413110-76fb452e-723b-4e74-b3ac-5952d54d663d.png)