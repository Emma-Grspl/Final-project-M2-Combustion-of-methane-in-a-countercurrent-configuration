# Simulation of Counter-Current Methane Combustion (M2 Project)

This project simulates methane–air combustion in a counter-current configuration. 
Methane and air are injected face-to-face in opposite directions, while nitrogen is injected 
from the upper and lower boundaries to dilute and stabilize the system. The computational 
domain measures 2 mm × 2 mm and features mixed boundary conditions:
- left wall: slip boundary  
- right wall: open boundary  
- top and bottom: no-slip, except injection zones  

The goal of this work is to solve the flow, species transport, and thermal chemistry to analyze 
the combustion behavior and identify key physical quantities such as deformation rates, 
diffusion zones, and peak temperature during the reaction.

## Objectives of the Project
### - Solve the flow field
We first solve the Navier–Stokes equations using a first-order fractional step method  
to obtain the velocity components:
- u(x, y) (horizontal velocity)  
- v(x, y) (vertical velocity)  

From the velocity field, we compute the deformation rate on the left wall:

$
a = \frac{dv}{dy}
$

This quantity is essential for understanding fluid–wall interactions, shear forces, and flow regime.

---

### **2. Transport of nitrogen**
Nitrogen behaves as a **spectator species** in the combustion process:
- It is not consumed  
- It does not affect temperature  
- It is simply convected and diffused  

Once the velocity field is known, we solve the transport equation for nitrogen to determine  
the **diffusive penetration zone** along the left boundary.  
This information is relevant for analyzing mixing, boundary layer development, and process optimization.

---

### **3. Full species transport + combustion temperature**
Finally, we simulate the transport of **all chemical species** participating in the reaction,  
and compute the associated **temperature evolution**.  
This allows us to estimate the **maximum flame temperature** during combustion.

For an ideal methane–air flame, the expected peak temperature ranges between **1700 and 2500 K**,  
which serves as a reference to validate the simulation.

---

## 📁 Repository Structure

