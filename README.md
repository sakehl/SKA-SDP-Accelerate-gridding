# SKA-SDP-Accelerate-gridding

Implementing the gridding algorithm from the [SKA Science Data
Processor][SKA-SDP] in [Accelerate], a domain specific language for
high performance computing embedded in [Haskell].



# Installation
* Make sure you have a version 1.10.0 of hdf5 installed. (The development library).
(Other versions might, work but I didn't test with them)
* Also alter the "extra-lib-dirs:" and "include-dirs" in package.yaml to point to the right directories for hdf5. For my ubuntu this is: "extra-lib-dirs: /usr/lib/x86_64-linux-gnu/hdf5/serial" and "include-dirs: /usr/include/hdf5/serial"
* When building and executing on ubuntu, I ran into the following linking [bug] for libffi. To fix it, I created the symbolic link libffi.so.7 pointing to libffi.so.6 in the systems lib directory
* The general accelerate instalation is used, with the CPU backend. The backend is build with LLVM 7.0, so install that one aswell. See [Accelerate] for more details.

# HDF

[bug]: https://github.com/commercialhaskell/stack/issues/4150
[SKA-SDP]:    https://github.com/SKA-ScienceDataProcessor
[Accelerate]: https://www.acceleratehs.org
[Haskell]:    https://www.haskell.org