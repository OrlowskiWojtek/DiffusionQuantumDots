# Diffusion Quantum Monte Carlo for 2D Quantum Dots

A C++ implementation of the Diffusion Quantum Monte Carlo (DQMC) method with Fixed-Node Approximation for calculating ground and excited states of quantum dots in anisotropic harmonic potentials.

## Purpose

This program enables high-accuracy quantum mechanical simulations of electron behavior in 2D quantum dot systems. It is particularly useful for:

- Calculating ground state energies of interacting electrons in quantum dots
- Approximating excited states using the Fixed-Node technique
- Studying electron density distributions in various confinement potentials
- Investigating electron-electron correlation effects in mesoscopic systems

The implementation uses importance sampling, population control, and the Fixed-Node Approximation to overcome the fermion sign problem, allowing for accurate simulations of fermionic systems.

## Features

- Jastrow-Slater trial wavefunctions for accurate nodal surfaces
- Harmonic oscillator basis functions for quantum dot systems
- OpenMP parallelization for improved performance
- Statistical analysis tools for error estimation
- Visualization capabilities for electron density distributions

## Documentation

For detailed information about the theoretical background, implementation details, and usage instructions, please visit the [comprehensive documentation](https://orlowskiwojtek.github.io/DiffusionQuantumDots/).
