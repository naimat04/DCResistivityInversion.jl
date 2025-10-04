# DC Resistivity Inversion

A high-performance Julia code for 1D DC resistivity inversion.

## 🚀 Features
- Forward modeling using analytical solutions
- Sensitivity analysis via brute-force method
- Regularized least-squares inversion
- Visualization of results

## 📈 Performance
- **Fast convergence**: 10 iterations to RMS ~1.28e-13
- **Robust**: Handles divergence with adaptive damping
- **Accurate**: Perfect data fit achieved

## 🎯 Quick Start
```julia
include("examples/run_inversion.jl")