# specsimfts

Implementation of "Simulation of stationary functional time series with given spectral density".



## Installation

```R
library("devtools")
install_github("tomasrubin/specsimfts")
```

## Running the demos

### FTS defined through its harmonic Karhunen-Loeve expansion

The simulation is the fastest and simplest if the FTS dynamics is defined directly through its harmonic Karhunen-Loeve, i.e. if the spectral density operator at each frequency is given by its spectral decomposition with its harmonic eigenvalues and harmonic eigenfunctions.

```R
library("specsimfts")
demo("demo_harmonicKL")
```

### FTS defined by its spectral density operators

If the harmonic expansion is not available but we define the FTS dynamics directly by its spectral density operator at each frequency, the spectral density operator needs to be discretised and decomposed by the SVD algorithm at each frequency. This method allows to simulate any dynamics defined by the spectral density operator but scales badly in terms of computational time as the discretisation resolution increases.

```R
library("specsimfts")
demo("demo_plain_spec_density")
```
    
