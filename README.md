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

From the velocity field, we compute the strain rate on the left wall: $ a = \frac{\partial v}{\partial y} $

This is essential to understand fluid–wall interactions and shear dynamics.

---

## 2. Nitrogen Transport  
Nitrogen behaves as a spectator species:

- It is not consumed  
- It does not affect temperature  
- It undergoes convection and diffusion only  

Using the computed velocity field, we solve the transport equation for \(Y_{N_2}\) to determine the **diffusive penetration zone** along the left boundary.

---

## 3. Full Species Transport + Combustion Temperature  
We simulate the transport and reaction of all key species:

- CH₄  
- O₂  
- CO₂  
- H₂O  

The temperature field is then computed from the reaction terms.  
For a methane–air flame, the expected peak temperature is **1700–2500 K**, providing a reference for validation.

---

# 📁 Project Structure
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
│
├── notebook/
│   ├── 01_velocity_field.ipynb
│   ├── 02_transport_N2.ipynb
│   └── 03_combustion.ipynb
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

---

# ▶️ How to Run the Project

## 1. Install dependencies

pip install -r requirements.txt

## 2. Run notebooks

jupyter notebook

## 3. Run a full simulation without notebooks

python examples/run_full_simulation.py
