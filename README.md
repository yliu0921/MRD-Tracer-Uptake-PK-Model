# MRD-Tracer-Uptake-PK-Model
This repository contains a pharmacokinetic (PK) model designed to predict the in vivo tracer uptake of molecular residual disease (MRD). It includes all relevant code and documentation used in the development of the model.

The project code consists of the following parts:

1. **`main.m`**：这是项目的主脚本，用于设置实验参数、加载数据以及调用所有子模块。它控制了整个实验流程，从数据输入、预处理、模型调用到结果的可视化输出。This is the main script for setting experimental parameters, loading data, and calling all submodules.
   
2. **`differential_eq.m`**：定义用于模拟MRD的Tracer PK过程的微分方程模型。方程组主要用于模拟放射性示踪剂在体内的分布和代谢情况。通过调整方程参数，可以灵活改变模型的行为，以适应不同的实验条件和假设。Define a differential equation model for the tracer PK processes. The equations are used to simulate the distribution and metabolism of radioactive tracer in the body. By adjusting the parameters of the equation, the behavior of the model can be changed flexibly to adapt to different experimental conditions and assumptions.
   
3. **`plot_results.m`**：用于对模拟的结果进行可视化，包括生成用于论文中的各种图表和结果展示。结果可能包括时间-活性曲线、MRD uptake定量以及不同实验条件下的比较分析。Used to visualize the results of the simulation, including the generation of various charts and results for use in the paper. The results include time-activity curves, MRD uptake quantification and comparative analysis under different experimental conditions.

4. **`parameters.m`**：包含项目中使用的所有全局参数设置，例如示踪剂的放射活性、注射剂量、PET扫描的时间间隔以及动物模型的生理参数。你可以在此文件中灵活调整各项实验参数。All global parameter settings used in the project were included, such as radioactivity of tracers, injection dose, time interval of PET scan, and physiological parameters of animal models. You can flexibly adjust the experimental parameters in this file.
