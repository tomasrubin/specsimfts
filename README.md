# Spectral domain simulation of functional time series 

This is the accompanying R package `specsimfts` for the paper: Rubin and Panaretos (2020). "Simulation of stationary functional time series with given spectral density".

The package contains all the methods introduced in the aformentioned paper including all the presented examples as demos (see below) that are easy to use and modify by functional time series (FTS) practitioners.



# Installation

```{r}
library("devtools")
install_github("tomasrubin/specsimfts")
```

## Running the demos

### FTS defined through eigendecomposition of its spectral density operators

The simulation is the fastest and simplest if the FTS dynamics is defined directly through the eigendecomposition of its spectral density operators. The simulation method then essentially mimics the Cramer-Karhunen-Loeve expansion.

```{r}
library("specsimfts")
demo("demo_harmonicKL")
```

### FTS defined by its spectral density operators

If the eigendecomposition of spectral density operators is not available but we still define the FTS dynamics directly by its spectral density operator at each frequency, the spectral density operator needs to be discretised and decomposed by the SVD algorithm at each frequency. This method allows to simulate any dynamics defined by the spectral density operator but scales badly in terms of computational time as the discretisation resolution increases.

```{r}
library("specsimfts")
demo("demo_plain_spec_density")
```
    
# Contact

Tomas Rubin

tomas.rubin@gmail.com

www.tomasrubin.com

[linkedin.com/in/tomas-rubin/](https://www.linkedin.com/in/tomas-rubin/)
