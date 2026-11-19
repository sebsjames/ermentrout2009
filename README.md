# Ermentrout 2009 reaction diffusion system

The Keller-Segal equations applied to whisker barrel cortex.

This implements the system described in *Subbarrel Patterns in Somatosensory Cortical Barrels Can Emerge from Local Dynamic Instabilities* by Ermentrout, Simons and Land:

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000537

![A screenshot of a reaction diffusion system which makes a cross-shaped pattern](https://github.com/sebsjames/ermentrout2009/blob/main/screenshots/celtic_cross.png?raw=true)


## Dependencies

This program uses C++ modules. At the time of writing, this means that
you have to compile with the clang compiler (version 20 or above).

If you are using Debian or Ubuntu, the following `apt` command should
install the mathplot dependencies.

```bash
sudo apt install build-essential cmake git ninja-build  \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev \
                 libglfw3-dev libfreetype-dev clang-20 clang-tools-20 \
                 librapidxml-dev libhdf5-dev
```

If you're building on a Mac, you can refer to mathplot's [Mac
README](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#installation-dependencies-for-mac)
for help. You only need to obtain and build
[glfw3](https://github.com/sebsjames/mathplot/blob/main/README.build.mac.md#glfw3);
OpenGL and Freetype should already be installed by default.

## Building

To build and run the example:

```bash
# Clone this example
git clone git@github.com:sebsjames/ermentrout2009

# Bring in the submodules - sebsjames/mathplot, sebsjames/maths , sebsjames/maths_config and nlohmann/json
cd ermentrout2009
git submodule init
git submodule update

# Build it in a 'build' directory
mkdir build
cd build
CC=clang-20 CXX=clang++-20 cmake .. -GNinja
ninja
cd ..

# Run an example
./build/erm configs/erm.json

# Reproduce the examples shown in the paper figures
./build/erm configs/erm_2009_fig5Ai.json
./build/erm configs/erm_2009_fig5Aii.json
./build/erm configs/erm_2009_fig5Aiii.json
./build/erm configs/erm_2009_fig5Aiv.json
```

## Reproduced figures

This program reproduces some of Fig 5. A i-iv. I don't get quite the same results for Aii or Aiii (the 'Mercedes' and the 'baseball').

### Ai 'Coffeebean' Dn = 200, Dc = 35

![A screenshot of my reproduction of Fig 5 A i](https://github.com/sebsjames/ermentrout2009/blob/main/screenshots/fig5Ai.png?raw=true)

### Aii 'Mercedes' Dn = 100, Dc = 18

![A screenshot of my reproduction of Fig 5 A ii](https://github.com/sebsjames/ermentrout2009/blob/main/screenshots/fig5Aii.png?raw=true)

### Aiii 'Baseball' Dn = 60, Dc = 14

![A screenshot of my reproduction of Fig 5 A iii](https://github.com/sebsjames/ermentrout2009/blob/main/screenshots/fig5Aiii.png?raw=true)

### Aiv 'Bullseye' Dn = 8, Dc = 2

![A screenshot of my reproduction of Fig 5 A iv](https://github.com/sebsjames/ermentrout2009/blob/main/screenshots/fig5Aiv.png?raw=true)
