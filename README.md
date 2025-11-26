# Simulation of Counter-Current Methane Combustion  
### Master 2 Project — Plasma and fusion physics - Ecole Polytechnique

This project simulates methane–air combustion in a counterflow configuration. Methane is injected from the left boundary, air from the right boundary, and nitrogen is injected from the upper and lower boundaries to dilute and stabilize the system. The computational domain measures 2 mm × 2 mm and includes mixed boundary conditions:
- Left wall : Slip boundary
- Right wall : Open boundary
- Top/bottom wall : No slip boundary (except injections)

The goals of this work are to solve the flow field, species transport, and chemical reaction to analyze the combustion behavior.

## Project Objectives
### Solve the Flow Field  

We solve the 2D Navier–Stokes equations using a first-order fractional-step method to compute:

- u(x, y) — horizontal velocity  
- v(x, y) — vertical velocity  

From the velocity field, we compute the strain rate on the left wall: a = dv/dy. This is essential to understand fluid–wall interactions and shear dynamics.

### Nitrogen Transport  
Nitrogen behaves as a spectator species:
- It is not consumed  
- It does not affect temperature  
- It undergoes convection and diffusion only  

Using the computed velocity field, we solve the transport equation for YN2 to determine the diffusive penetration zone along the left boundary.

### Full Species Transport + Combustion Temperature  
We simulate the transport and reaction of all key species:
- CH4  
- O2 
- CO2  
- H20  
The temperature field is then computed from the reaction terms. For a methane–air flame, the expected peak temperature is 1700–2500 K, providing a reference for validation.

### Final_Project_Combusion_Methane.ipynb

All the details of this work with the single-block code are in this notebook. If you want more information on the physics involved, I encourage you to read it. In particular, to understand all the equations and notations used throughout this project

## Project Structure
```bash
combustion-counterflow/
│
├── src/
│   └── counterflow/
│       ├── __init__.py
│       ├── config.py
│       ├── grid.py
│       ├── navier_stokes.py
│       ├── transport.py
│       └── combustion.py
│       └── analysis.py
│
├── notebook/
│   ├── 01_velocity_field.ipynb
│   ├── 02_transport_N2.ipynb
│   └── 03_combustion.ipynb
│   └── 04_analysis.ipynb
│   └── raw notebook/
│       └── Final_Project_Combustion_Methane
│
├── examples/
│   └── run_full_simulation.py
│
├── tests/
│   ├── test_navier_stokes.py
│   ├── test_transport.py
│   └── test_combustion.py
│
├── figures/
│
├── requirements.txt
└── LICENSE

```

## How to Run the Project

### Install dependencies

pip install -r requirements.txt

### Run notebooks

jupyter notebook

### Run a full simulation without notebooks

python examples/run_full_simulation.py
