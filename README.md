# Solver for the Orr-Sommerfeld equation

Solver of the Orr Sommerfeld equation, both for the  temporal and spatial stability analysis. The code follows the work of: 

_Kirchner, N.P. (2000), Computational aspects of the spectral Galerkin FEM for the Orr–Sommerfeld equation. Int. J. Numer. Meth. Fluids, 32: 105-121. [https://doi.org/10.1002](https://doi.org/10.1002/(SICI)1097-0363(20000115)32:1<105::AID-FLD938>3.0.CO;2-X)_

The paper can be found in `docs/spectralElementMethodOrrSommerfeld_0.pdf`. The Galerkin method has the advantage of avoiding the appearance of spurious eigenvalues, which is a common issue in collocation-based spectral codes.

## Requirements

- Python 3.7 or later
- cmake
- gcc
- Eigen library

## Installation

Copy this repository to your local machine:

```bash
git clone git@github.com:victorballester7/orrSommerfeldSolver.git
cd orrSommerfeldSolver
make
```

## Usage

All the parameters to change are in the `config/input.toml` file. The comments should be self-explanatory. Once the user has selected the parameters, the user can run the code with:

```bash
make run
```


