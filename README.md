# PyCantPhi

PyCantPhi is a Python package for calculating anthropogenic carbon (Cant) in the North Atlantic Ocean using the φ CT° (Phi) method, and cStar and Troca methods. The cStar and Troca methods can be applyied globally for general purpose.

## Features

- Implementation of the φ CT° method for anthropogenic carbon calculation
- Integration with Mauna Loa CO2 data
- Support for water mass analysis using Optimum MultiParameter (OMP)
- Preformed alkalinity calculations
- Mediterranean water influence calculations
- CO2 system calculations using PyCO2SYS

## Installation

You can install PyCantPhi using pip:

```bash
pip install pycantphi
```

## Quick Start

```python
import numpy as np
import pandas as pd
import xarray as xr
from pycantphi import cantphi, CantCalculator


# Define dimensions for dataset
loc = ["loc1", "loc2", "loc3"]
pressure = [0, 250, 2000]  # in dbar
time = pd.date_range(start='2020', end='2023', freq='Y')
lon = [-10, -20, -40]
lat = [10, 20, 65]

# Create coordinate data for xarray
coords = {
    "location": loc,
    "pressure": pressure,
    "time": time,
    "longitude": ("location", lon),
    "latitude": ("location", lat), 
    "year": ("time", time.year),   
}

# Note that it's also possible to have no location variable (here, we're creating a three-dimensional dataset) and rely solely on longitude or latitude. The other variable (longitude or latitude) is always required as a function of the first. 

# Generate synthetic data arrays 
data_vars = {
    "theta": (["location", "pressure", "time"], np.random.uniform(1, 15, (3, 3, 3))),  # in Celsius
    "salinity": (["location", "pressure", "time"], np.random.uniform(30, 37, (3, 3, 3))),    # in PSU
    "alkalinity": (["location", "pressure", "time"], np.random.uniform(2200, 2400, (3, 3, 3))),  # in µmol/kg
    "oxygen": (["location", "pressure", "time"], np.random.uniform(310, 330, (3, 3, 3))),  # in µmol/kg
    "carbon": (["location", "pressure", "time"], np.random.uniform(1800, 2100, (3, 3, 3))), # in µmol/kg
    "phosphate": (["location", "pressure", "time"], np.random.uniform(0, 3, (3, 3, 3))),   # in µmol/kg
    "nitrate": (["location", "pressure", "time"], np.random.uniform(0, 40, (3, 3, 3))),    # in µmol/kg
    "silicate": (["location", "pressure", "time"], np.random.uniform(0, 150, (3, 3, 3)))   # in µmol/kg
}

# Create the xarray Dataset
ds = xr.Dataset(data_vars, coords=coords)

# Display the dataset
ds

# Initialize the calculator
calc = cantphi(ds)

# Process the dataset
results = calc.process()

# Access the calculated Cant values
cant_phi = results.cAntPhiCt0ML
cant_troca = results.cAntTroca
cant_cstar = results.cStar

# Access only the Cant calculator for the different methods. They can be run outside but required specific parameters. 

CantCalculator.calculate_cstar(ct, aou, ct_preformed, delta_ca_corrected)
CantCalculator.calculate_ctroca(ct, at, oxygen, theta)
CantCalculator.calculate_cphi(calc.ds)

```

## Required Data Format

Your input dataset should include the following variables:
- longitude
- latitude
- pressure
- theta (potential temperature)
- salinity
- oxygen
- silicate
- nitrate
- phosphate
- carbon
- alkalinity
- year

## Documentation

The documentation is in progress. For detailed documentation, examples, and API reference, visit our [documentation page](https://pyphi.readthedocs.io) once finished.
```

## Contributing

Contributions are welcome! I'm working on generating pytests before taking new contributions into account. If you encounter a problem with the package, please raise an issue.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite both the original CantPhiCt0 paper and this implementation in python:

For the method:

```
@article{vázquez2009anthropogenic,
  title={Anthropogenic carbon distributions in the Atlantic Ocean: data-based estimates from the Arctic to the Antarctic},
  author={V{\'a}zquez-Rodr{\'\i}guez, M and Touratier, F and Lo Monaco, C and Waugh, DW and Padin, XA and Bellerby, RGJ and Goyet, C and Metzl, N and R{\'\i}os, AF and P{\'e}rez, FF},
  journal={Biogeosciences},
  volume={6},
  number={3},
  pages={439--451},
  year={2009},
  doi={10.5194/bg-6-439-2009},
  publisher={Copernicus GmbH}
}
```

For the software:

```
@software{pyphi2024,
  title = {PyCantPhi: A Python implementation of the φ CT° method for anthropogenic carbon calculation},
  author = {Bajon R.},
  year = {2024},
  version = {0.0.0},
  url = {https://github.com/RaphaelBajon/pyphi}
}
```

## Acknowledgments

The implementation is based on the φ CT° method described in Vázquez-Rodríguez et al. (2009).

