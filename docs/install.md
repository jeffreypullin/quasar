# Installation

The quasar source code is hosted on [GitHub](https://github.com/jeffreypullin/quasar)

# Installation

To compile quasar from source run:

```sh
git clone https://github.com/jeffreypullin/quasar.git
cd quasar
mkdir -p build
cd build
cmake ..
make
```

This process generates the binary in the `build/` subdirectory of quasar. To verify the installation has completed succesfully, run

```
./quasar --version
```
