# Greedy Align

This repository contains experimental code for the greedy alignment algorithm.

## Installation

This tool is built with [SeqAn3](https://docs.seqan.de/seqan/3-master-user/index.html). To properly build the package, you need to have GCC >= 11.3, G++ and CMake installed.

```bash
git clone https://github.com/GZHoffie/bucket-map.git
cd bucket-map
mkdir build
cd build

cmake ../
```

This will automatically download the dependencies. Next, we can build the targets. For example, for testing using the main file [`test.cpp`](./test.cpp), we can build using 

```bash
cmake --build . --target greedy_aligner_test
```

And for benchmarking the speed of aligner, we can build using the main file [`main.cpp`](./main.cpp) (currently buggy).

```bash
cmake --build . --target greedy_aligner
```

