# HCRefrigeratorModel
This repository provides a structured and computationally efficient tool for researchers and engineers exploring hydrocarbon refrigerants and their impact on refrigeration system performance. The repository contains the simulation code for the paper "Hydrocarbon mixtures as Alternative Refrigerants in Household Refrigerators: Dynamic Modeling, Experimental Validation, and Parametric Study" published in Scientific Reports journal with doi code: xxxx.

This repository provides a simulation framework for studying refrigerators that use two-component hydrocarbon refrigerant mixtures. It combines a fast, regression-based thermodynamic model with detailed component-level calculations, making it possible to analyze system performance efficiently and reliably. The project begins with a base case using pure R600a refrigerant and then extends to binary mixtures such as R290:R600a, R290:R600, R290:R1270, R600a:R600, R600a:R1270, and R600:R1270. Within the code, refrigerant properties are estimated using regression for speed, and the refrigeration cycle is modeled through its major components.

# Thermodynamic model
This project focuses on simulating a refrigerator that operates with two-component hydrocarbon refrigerant mixtures. The mixtures under study include:
- R290:R600a
- R290:R600
- R290:R1270
- R600a:R600
- R600a:R1270
- R600:R1270

A fast regression-based thermodynamic model is employed for hydrocarbon property calculations. This approach significantly reduces computational time, making the model suitable for optimization tasks where repeated calculations are required.
# Refrigerator simulation
The simulation is carried out in two main phases:
1. Basecase Simulation with Pure Refrigerant
   - The system is first modeled using pure R600a refrigerant to establish a reference case.
2. Case Study with Binary Mixtures
   - A detailed analysis is then performed for the selected binary hydrocarbon mixtures listed above.
# Structure of the Code
The code follows a systematic approach:
1. Thermodynamic Data Retrieval and Regression
   - Refrigerant thermodynamic properties are obtained and approximated using regression techniques for computational efficiency.
2. Component-Level Calculations
   - The refrigeration cycle is broken down into major components, with calculations performed in sequence for:
     - Evaporator: heat absorption process
     - Cabinet model: thermal load representation of the refrigerator cabinet
     - Compressor model: work input and compression process
3. Technical Performance Calculations
   - Key performance metrics (e.g. compressor work) are computed.
4. Numerical Solution of ODEs
   -Taylor series expansion is applied to solve the ordinary differential equations (ODEs) arising from the system’s dynamic or component-level modeling.
