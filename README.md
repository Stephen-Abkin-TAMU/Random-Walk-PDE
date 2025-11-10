# Random-Walk-PDE
A Python framework for simulating PDEs (Heat, FitzHugh-Nagumo, Burgers) using gradient random walks and deterministic schemes.

The project demonstrates how random walk approaches can replicate diffusion and reaction–diffusion dynamics while also providing a comparison with traditional grid-based solvers.

  ## Features
  - Supports multiple PDEs: Heat equation, FitzHugh–Nagumo model, and Burgers’ equation
  - Gradient Random Walk method for particle-based diffusion without fixed grids
  - Deterministic finite-difference solver for Burgers’ equation
  - Configurable parameters: domain type, boundary conditions, diffusivity, time step, total time, number of points
  - Generates plots of simulation results, saved automatically in the `output/` folder

  ## Applications
  - Heat flow and diffusion modeling
  - Nerve signal propagation
  - Fluid dynamics and nonlinear PDEs
  - Exploring trade-offs between stochastic and deterministic solvers
