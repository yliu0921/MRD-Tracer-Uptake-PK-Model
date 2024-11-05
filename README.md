# MRD-Tracer-Uptake-PK-Model
This repository contains a pharmacokinetic (PK) model designed to predict the in vivo tracer uptake of molecular residual disease (MRD). It includes all relevant code and documentation used in the development of the model.

The project code consists of the following parts:

1. **`FigureX.m`**：These are the main scripts for setting experimental parameters, loading data, and calling all submodules. All global parameter settings used in the project were included, such as radioactivity of tracers, injection dose, time interval of PET scan, and physiological parameters of animal models. You can flexibly adjust the experimental parameters in this file. Used to visualize the results of the simulation, generation of each figures used in the paper. The results include time-activity curves, MRD uptake quantification and comparative analysis under different experimental conditions.
   
2. **`ODEs.m`**：Define the ordinary differential equations for the tracer PK processes in tumour microenviroment. By adjusting the parameters of the equation, the behavior of the model can be changed flexibly to adapt to different experimental conditions and assumptions.
   
3. **`Inputs.m`**：The initial conditions (y0) and parameter matrix (P) for the ODE system defined above were calculated by passing molecular weight(MW), number of receptors on tumour cell surface(NR), binding constant(kon), dissociation constant(koff) and endocytosis rate(kendo). One can also modify the system behavior by passing different parameter values through this file to mimic different experimental conditions. The interstitium void fraction 'SchmidtVoid.m' and permeability 'SchmidtPerm.m' calculation method from Schmidt's study are embedded in this function to further characterize the tumor microenvironment.

4. **`Cal_Rk.m`**：Calculate the microlesion size based on the characteristics of the tumour cells, such as cell diameter, cell number, and the tumor interstitium fluid fraction.

5. **`Cal_kcl.m`**：Calculate the plasma clearance rate of the tracers based on its molecular weight, with the conversion relationship derived from a nonlinear cubic fit of empirical values extracted from multiple studies published previously.

### Example Process

1. **Set global parameters in `FigureX.m`**;
2. **Call the `Inputs` function in `main.m`**: Call the `Inputs` function in `FigureX.m` to obtain the ODE initial conditions and parameter matrix based on the global parameters set;
3. **Run the ODE solver**: Solve the differential equations using the returned parameter matrix.

---
