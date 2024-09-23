# MRD-Tracer-Uptake-PK-Model
This repository contains a pharmacokinetic (PK) model designed to predict the in vivo tracer uptake of molecular residual disease (MRD). It includes all relevant code and documentation used in the development of the model.

The project code consists of the following parts:

1. **`main.m`**：这是项目的主脚本，用于设置实验参数、加载数据以及调用所有子模块。它控制了整个实验流程，从数据输入、预处理、模型调用到结果的可视化输出。包括生成用于论文中的各种图表和结果展示。结果包括时间-活性曲线、MRD uptake定量以及不同实验条件下的比较分析。This is the main script for setting experimental parameters, loading data, and calling all submodules.Used to visualize the results of the simulation, including the generation of various charts and results for use in the paper. The results include time-activity curves, MRD uptake quantification and comparative analysis under different experimental conditions.
   
2. **`differential_eq.m`**：定义用于模拟MRD的Tracer PK过程的微分方程模型。方程组主要用于模拟放射性示踪剂在体内的分布和代谢情况。通过调整方程参数，可以灵活改变模型的行为，以适应不同的实验条件和假设。Define a differential equation model for the tracer PK processes. The equations are used to simulate the distribution and metabolism of radioactive tracer in the body. By adjusting the parameters of the equation, the behavior of the model can be changed flexibly to adapt to different experimental conditions and assumptions.
   
3. **`Inputs.m`**：传递了 MW（分子量）、NR（肿瘤细胞上的受体密度）、结合常数、解离常数、内吞速率等参数，来计算 ODE 系统所需的初始条件 (y0) 和参数矩阵 (p)。
您也可以根据不同的实验条件通过`parameters.m`传入不同的参数值来改变系统行为。The initial conditions (y0) and parameter matrix (P) for the ODE system were calculated by passing MW, NR, binding constant, dissociation constant and endocytosis rate.You can also modify the system behavior by passing different parameter values through `parameters.m` based on different experimental conditions.

4. **`parameters.m`**：包含项目中使用的所有全局参数设置，例如示踪剂的放射活性、注射剂量、PET扫描的时间间隔以及动物模型的生理参数。你可以在此文件中灵活调整各项实验参数。All global parameter settings used in the project were included, such as radioactivity of tracers, injection dose, time interval of PET scan, and physiological parameters of animal models. You can flexibly adjust the experimental parameters in this file.

### 1 Example Process

1. **Set global parameters in `parameters.m`**: Set and manage the key experimental parameters centrally through `parameters.m`.
2. **Call the `Inputs` function in `main.m`**: Call the `Inputs` function in `main.m` to obtain the ODE initial conditions and parameter matrix based on the global parameters set.
3. **Run the ODE solver**: Solve the differential equations using the returned `y0` and `p`.


### 2 Parameters Explanation
This section provides detailed explanations of the parameters defined in the parameters.m file. These parameters control various aspects of the system behavior, including tracer dose, binding constants, cell properties, and numerical simulation settings. By adjusting these parameters, you can modify the conditions under which the system is simulated.

#### Tracer and Cell-Related Parameters

- **`params.L`**: The tracer dose, used as the initial concentration of the tracer in the simulation.
- **`params.LIDcircV`**: The concentration of the tracer in the body, representing molecules per unit volume.
- **`params.NR`**: The receptor density on tumor cells, indicating the number of receptors per cell.

#### PET Scan Settings

- **`params.Vb`**: The volume of distribution of the drug in the body (in liters), used to model how the drug spreads.
- **`params.scan_time_interval`**: The time interval between PET scans, specified in hours. This controls how frequently the system is scanned.

#### Binding and Dissociation Constants

- **`params.kon_R`**: The on-rate binding constant for the interaction between the receptor and drug, measured in 1/M/s. This parameter controls how quickly the drug binds to the receptor.
- **`params.koff_R`**: The off-rate dissociation constant, representing how quickly the drug-receptor complex dissociates.
- **`params.kendo_R`**: The endocytosis rate constant, specifying the rate at which the receptor-drug complex is internalized by cells.
- **`params.Lesc`**: The initial concentration of the escape factor, representing the ability of the drug to diffuse or escape from the tumor environment.

#### ODE Solver Settings

- **`params.tspan`**: The time span for the ODE solver, indicating the total duration of the simulation (in seconds). This defines how long the system is simulated.
- **`params.tplot`**: Time points used for plotting the results, allowing for detailed visualization of the system behavior at various time intervals.
- **`params.options`**: ODE solver accuracy options, including relative and absolute tolerances to control the precision of the numerical solution.

#### Numerical Simulation Iteration Parameters

- **`params.elem`**: The number of elements used in numerical iteration, representing how many steps or data points are considered in the simulation.
- **`params.koff_R_array`**: An array of dissociation constants used to explore the behavior of the system under different dissociation rates.
- **`params.MW_array`**: An array of molecular weights used to simulate the system for different drug sizes.
- **`params.cellNum_array`**: An array of cell numbers, allowing for simulations with varying tumor cell densities.

#### Specific Experimental Parameters

- **`params.cellNum`**: The number of tumor cells in the system, influencing the total number of receptors present in the simulation.

#### Physical Parameters for the System

- **`params.Mol_R_factor`**: A factor used to calculate the molecular radius based on molecular weight.
- **`params.Rcap`**: The capillary radius, which is important for modeling drug diffusion and transport within the tumor microenvironment.
- **`params.Rcell`**: The radius of a tumor cell, which impacts how the drug interacts with individual cells.
- **`params.Vi`**: The intracellular volume factor, used to calculate the intracellular distribution of the drug.
- **`params.Rk`**: The Krogh cylinder radius, which defines the geometry of the tissue region being modeled.
- **`params.P_factor`**: The permeability constant, controlling how easily the drug can permeate through the tumor tissue.

---
