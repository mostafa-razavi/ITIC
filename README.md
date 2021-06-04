# ITIC
The ITIC method to calculate VLE data of pure fluids from MD or MC simulation data 


# Project Title

Simple overview of use/purpose.

## Description

The isothermal-isochoric integration (ITIC) method is a viable method for vapor-liquid coexistence calculation by molecular simulation. The ITIC method proves to be much more effective compared to GEMC and GCMC methods for vapor-liquid coexistence calculations at reduced temperatures of 0.45–0.6, which are important for practical applications. Furthermore, the ITIC method lends itself to application with molecular dynamics (MD) as well as MC, advancing the prospect of simulation results that are quantitatively consistent across software platforms.

### Dependencies

* gfortran

### Installing

* Download the repository on your local machine
* Compile the ITIC executable file by running bash compile.sh in src/ directory. You will need gfortran installed on your computer

### Executing program

* Use the ITIC_selector.xlsx spreadsheet to generate an ITIC input file for your molecule
* Run NVT molecular simulation (MD or MC) to obtain Z and Ures data at the ITIC points determined in the previous step
* Process simlation data into a data file with a format similar to example files (see C4/C4_MiPPE.zures for example)
* Modify the first four variables in the ITIC.sh script (see C4/ITIC.sh for example)
* Run bash ITIC.sh


If you have a question or comment, please feel free to contact S. Mostafa Razavi at sr87@uakron.edu


## Cite the ITIC paper
Razavi, S. Mostafa, Richard A. Messerly, and J. Richard Elliott. "Coexistence calculation using the isothermal-isochoric integration method." Fluid Phase Equilibria 501 (2019): 112236.
https://doi.org/10.1016/j.fluid.2019.06.026

Bibtex citation record:

@article{razavi2019coexistence,
  title={Coexistence calculation using the isothermal-isochoric integration method},
  author={Razavi, S Mostafa and Messerly, Richard A and Elliott, J Richard},
  journal={Fluid Phase Equilibria},
  volume={501},
  pages={112236},
  year={2019},
  publisher={Elsevier}
}

