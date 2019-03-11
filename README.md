# SKA-SDP-Accelerate-gridding

Implementing the gridding algorithm from the [SKA Science Data
Processor][SKA-SDP] in [Accelerate], a domain specific language for
high performance computing embedded in [Haskell].

[SKA-SDP]:    https://github.com/SKA-ScienceDataProcessor
[Accelerate]: https://www.acceleratehs.org
[Haskell]:    https://www.haskell.org

# Installation
* Make sure you have a version 1.10.0 of hdf5 installed. (The development library).
(Other versions might, work but I didn't test with them)
* Also alter the "extra-lib-dirs:" and "include-dirs" in package.yaml to point to the right directories for hdf5. For my ubuntu this is: "extra-lib-dirs: /usr/lib/x86_64-linux-gnu/hdf5/serial" and "include-dirs: /usr/include/hdf5/serial"
