# fts-spectral-simulation
Toolbox for spectral domain simulation of functional time series samples

print(paste("Nuclear norm relative error:", sum(svd(r, nu=0, nv=0)$d) / sum(diag(covlag0))))