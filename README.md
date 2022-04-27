# ELMA-BraCat
ELMA - Extensions of the LeMonADE library
BraCat - Branched Catalysis related to the [paper](https://pubs.acs.org/doi/abs/10.1021/jacs.9b06785)

## Installation

* Clone and Install LeMonADE v2.2.2, see [LeMonADE v2.2.2](https://github.com/LeMonADE-project/LeMonADE/tree/v2.2.2)
* Install cmake (minimum version 2.8)
* Just do for standard compilation of the ELMA-BraCat project:
````sh
    # generates the projects
    mkdir build
    cd build
    cmake -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/ ..
    make
````
or

````sh
    cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/
    cmake --build build/
````
* The executables can be found in the build/projects folder.

## Environment
* tested with g++ version GNU 11.1.1
* [LeMonADE v2.2.2](https://github.com/LeMonADE-project/LeMonADE/tree/v2.2.2) or [LeMonADE v2.2.2 on Zenodo](https://doi.org/10.5281/zenodo.5061542)
* cmake version 3.20.4

## Run  

### Create branched structures with Chain Walking (CW) approach

* Creates an excluded volume CW structure with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box:
````sh
    ./CreatorChainWalkingTertiaryBondWalking -f out.bfm -n 500 -p 0.5 -t 1.0 -b 256
````
* Evaluates the creation properties (distance between reaction events) of an excluded volume CW structure with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box using 10 samples:
````sh
    ./CreatorChainWalkingTertiaryBondWalking_CreationProperties -f out.bfm -n 500 -p 0.5 -b 256 -s 10
````
* Evaluates the topological properties (number vertices from topological center) of an excluded volume CW structure with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box using 10 samples:
````sh
    ./CreatorChainWalkingTertiaryBondWalking_TopologicalProperties -f out.bfm -n 500 -p 0.5 -b 256 -s 10
````
* Evaluates the eigenvalues of the Rouse matrix of an excluded volume CW structure with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box using 10 samples:
````sh
    ./CreatorChainWalkingTertiaryBondWalking_RouseMatrix -f out.bfm -n 500 -p 0.5 -b 256 -s 10
````
* Simulation of a BFM structure under excluded volume condition for 1000MCS with dump every 100MCS:
````sh
    ./SimpleSimulatorStandard -i in.bfm -o out.bfm -m 1000 -s 100
````
* Creates an excluded volume Barabasi-Albert structure with restricted node functionality $f <= 3$ with 500 monomers, $m_0=2$, and $m=1$ in $L=256$ box:
````sh
    ./CreatorBarabasiAlbertRestrictedF3_Creator -f out.bfm -n 500 -b 256
````
* Evaluates the eigenvalues of the Rouse matrix of an excluded volume Barabasi-Albert structure with restricted node functionality $f <= 3$ with 500 monomers, $m_0=2$, and $m=1$ in $L=256$ box using 10 samples:
````sh
    ./CreatorBarabasiAlbertRestrictedF3_RouseMatrix -f out.bfm -n 500 -b 256 -s 10
````
* Evaluates the creation properties (distance between reaction events) of ideal BGRW structure (no excluded volume) with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box using 10 samples:
````sh
    ./CreatorBGRWTertiaryBondWalking_CreationProperties -f out.bfm -n 500 -p 0.5 -b 256 -s 10
````
* Evaluates the eigenvalues of the Rouse matrix of ideal BGRW structure (no excluded volume) with 500 monomers, walking rate $w=1/p=1/0.5=2$ in $L=256$ box using 10 samples:
````sh
    ./CreatorBGRWTertiaryBondWalking_RouseMatrix -f out.bfm -n 500 -p 0.5 -b 256 -s 10
````
* Creates an ideal Barabasi-Albert structure (no excluded volume) with 500 monomers, $m_0=2$, and $m=1$ in $L=256$ box:
````sh
    ./CreatorBarabasiAlbert_Creator -f out.bfm -n 500 -b 256
````
* Evaluates the eigenvalues of the Rouse matrix of an ideal Barabasi-Albert structure (no excluded volume) with 500 monomers, $m_0=2$, and $m=1$ in $L=256$ box using 10 samples:
````sh
    ./CreatorBarabasiAlbert_RouseMatrix -f out.bfm -n 500 -b 256 -s 10
````
* Evaluates the time average of the radius of gyration $Rg2$ of any structure:
````sh
    ./ChainWalking_Analyzer_RG2 -f out.bfm
````
* Evaluates the eigenvalues of the Rouse matrix of an excluded volume Slow Growth structure with restricted node functionality $f <= 3$ with 500 monomers in $L=256$ box using 10 samples:
````sh
    ./CreatorSlowGrowth_RouseMatrix -f out.bfm -n 500 -b 256 -s 10
````
* Creates a perfect dendrimer structure under excluded volume condition with functionality $f <= 3$, generation $G=4$, and spacer length $S=2$ in $L=256$ box:
````sh
    ./CreateDendrimerGeneral_FGS -o out.bfm -f 3 -g 4 -s 2 -x 256 -y 256 -z 256
````
* Evaluates the eigenvalues of the Rouse matrix of an perfect "hyperstar" structure (no excluded volume - $2^(-f)$ node distribution) with generation $G=7$ in $L=256$ box:
````sh
    ./CreatorHyperstar_RouseMatrix -f out.bfm -g 7 -b 512
````

## License

See the LICENSE in the root directory.

## How to cite

The development of the code is mainly funded by academic research grants.   
If you find the source code useful please cite the related paper:

[1] A. Jurjiu, R. Dockhorn, O. Mironova, J.-U. Sommer,  
    "Two universality classes for random hyperbranched polymers", [Soft Matter 2014, 10, 4935-4946](https://doi.org/10.1039/C4SM00711E)  
[2] R. Dockhorn, L. Plüschke, M. Geisler, J. Zessin, P. Lindner, R. Mundil, J. Merna, J.-U. Sommer, A. Lederer,  
    "Polyolefins Formed by Chain Walking Catalysis-A Matter of Branching Density Only?", [J. Am. Chem. Soc. 2019, 141 (31), 15586-15596](https://doi.org/10.1021/jacs.9b06785)  
[3] R. Dockhorn, J.-U. Sommer,  
    "Theory of Chain Walking Catalysis: From Disordered Dendrimers to Dendritic Bottle-Brushes", 2022
