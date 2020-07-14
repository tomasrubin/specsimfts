# Spectral domain simulation of functional time series 

This is the accompanying R package `specsimfts` for the paper:

**Rubin and Panaretos (2020). "Simulation of stationary functional time series with given spectral density".**

The package contains all the methods introduced in the aformentioned paper including all the presented examples as demos (see below) that are easy to use and modify by functional time series (FTS) practitioners.


## Methods description

The packages includes four simulation approaches (see below). Each approach includes a function to simulate a FTS sample and another one to calculate the theoretical lagged autocovariance operator. The numerical calculation of these operators is performed by integrating the spectral density operators through the inverse formula. While the simulation methods are generally fast, the inverse formula integration can be slower.

* Simulation using the Cramer-Karhunen-Loeve expansion:
    + functions: `CKL_simulate`, `CKL_covlagh_operator`,
    + demo (Example 4.1 in the paper): `demo_CKL`.
* Simulation with given spectral density operators (with no special structure),
    + functions: `spec_density_simulate`, `spec_density_covlagh_operator`,
    + demo (Example 4.1 in the paper): `demo_plain_spec_density`.
* Simulation of filtered white noise:
    + functions: `filter_simulate`, `filter_covlagh_operator`,
    + demos (Example 4.2 in the paper): `demo_FARFIMA_as_filter`, `demo_FARFIMA_as_filter_SVD`,
    + additional demos for custom-defined filter: `demo_custom_filter`, `demo_custom_filter_SVD`.
* Simulation of general FAR(FI)MA processes:
    + functions: `FARFIMA_test_stationarity`, `FARFIMA_simulate`, `FARFIMA_covlagh_operator`,
    + demos (Example 4.3 in the paper): `demo_FARMA`, `demo_FARMA_SVD`,
    + additional demo: `demo_FARFIMA_spec_vs_hybrid`.
    
See below on guidence on how to run the demos and to read their description.

## Installation

```{r}
library("devtools")
install_github("tomasrubin/specsimfts")
```

## Running the demos

### FTS specified through eigendecomposition of its spectral density operators

The simulation is the fastest and simplest if the FTS dynamics is defined directly through the eigendecomposition of its spectral density operators. The simulation method then essentially mimics the Cramer-Karhunen-Loeve expansion.

```{r}
library("specsimfts")
demo("demo_CKL")
```

### FTS specified by its spectral density kernels

If the eigendecomposition of spectral density operators is not available but we still define the FTS dynamics directly by its spectral density operator at each frequency, the spectral density operator needs to be discretised and decomposed by the SVD algorithm at each frequency. This method allows to simulate any dynamics defined by the spectral density operator but scales badly in terms of computational time as the discretisation resolution increases.

```{r}
library("specsimfts")
demo("demo_plain_spec_density")
```
###  FTS specified as filtered white noise

The FARFIMA(1,0.2,0) process scrutinised in Example 4.2 can be written as a filtered white noise process while having a direct formula for the frequency response function because of the special structure of the autoregressive operator. The demo file `demo_FARFIMA_as_filter.R` contains this implementation together with benefiting from the analytic eigendecomposition of the innovation white noise processes.

```{r}
library("specsimfts")
demo("demo_FARFIMA_as_filter")
```

If the innovation white noise processes does not admit analytically known eigendecomposition it can be calculated by the SVD algorithm, inspect the demo below.

```{r}
library("specsimfts")
demo("demo_FARFIMA_as_filter_SVD")
```

As an additional demo file we provide with an extra example where the frequency response function is artificially defined using operations including rank-one-tensors, kernel integral operators, antiderivatives, and identity operators. See the comments in the demo file. The white noise admits Brownian motion covariance with known analytic eigendecomposition

```{r}
library("specsimfts")
demo("demo_custom_filter")
```

or which can be numerically calculated by the SVD algorithm

```{r}
library("specsimfts")
demo("demo_custom_filter_SVD")
```

### FAR(FI)MA processes with general autoregressive and moving average operators (integral operators specified by their kernels)

The FARMA(4,3) process, being a special case of the FARFIMA(p,d,q) process with d=0, scrutinized in Example 4.3 in the paper is included as a demo file `demo_FARMA`. The code allows to include an arbitrary number of autoregressive and moving average operators defined as integral operators with given kernels, and to set the fractional integration parameter d in the interval (-0.5, 0.5).

```{r}
library("specsimfts")
demo("demo_FARMA")
```

While the above example uses the explicite finite-rank specification of the innovation white noise, the FAR(FI)MA processes can be also simulated with an arbitrary innovation noise whose eigendecomposition is numerically calculated by the SVD algorithm.

```{r}
library("specsimfts")
demo("demo_FARMA_SVD")
```

Besides the fully spectral method which can be a bit slower for FARFIMA processes with non-trivial autoregressive part, we have also proposed a hybrid simulation method. The demo file below shows how to run the two approaches and compares their running times.

```{r}
library("specsimfts")
demo("demo_FARFIMA_spec_vs_hybrid")
```

## Usage 

Individuals are free to use the codes for the purpose academic research, provided it is properly acknowledged and the paper *Rubin and Panaretos (2020) "Simulation of stationary functional time series with given spectral density"* is cited. For any other use, permission must first be arranged with the author. Unless otherwise specified, the author of the codes is Tomas Rubin (tomas.rubin@gmail.com). Please contact me if you find errors in the codes or if you have feature requests.
    
## Contact

Tomas Rubin

tomas.rubin@gmail.com

www.tomasrubin.com

[linkedin.com/in/tomas-rubin/](https://www.linkedin.com/in/tomas-rubin/)
