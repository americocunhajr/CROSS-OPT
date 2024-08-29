## Cross-Entropy Optimization for Truss Structures

**CROSS-OPT: Cross-Entropy Optimization for Truss Structures** is a Matlab package that implements a framework for size and shape structural optimization of truss systems. The package leverages the Cross-Entropy (CE) method for global search optimization and employs an augmented Lagrangian formulation to handle equality and inequality constraints. With some straightforward adaptations, the CROSS-OPT strategy can also be applied to other structural systems and optimization problems. 

<p align="center">
<img src="logo/CROSS-OPTframework.png" width="80%">
</p>

**CROSS-OPT** uses as optimization tool the package **CEopt - Cross-Entropy Optimizer**, which can be downloaded at <a href="https://ceopt.org" target="_blank">https://ceopt.org</a>.

### Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Usage](#usage)
- [Documentation](#documentation)
- [Reproducibility](#reproducibility)
- [Authors](#authors)
- [Citing CROSS-OPT](#citing-cross-opt)
- [License](#license)
- [Institutional support](#institutional-support)
- [Funding](#funding)

### Overview
**CROSS-OPT** was developed to solve nonconvex optimization problems involving truss structures. The underlying results are reported in the following publication:
More details are in the following paper:
- *M. V. Issa, A. Pereira and A. Cunha Jr, The cross-entropy method for nonconvex structural optimization, 2024 (under review)*

Preprint available <a href="https://hal.archives-ouvertes.fr/xxx" target="_blank">here</a>.

### Features
- 

### Usage
To get started with **CROSS-OPT**, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/americocunhajr/CROSS-OPT.git
   ```
2. Navigate to the code directory:
   ```bash
   cd CROSS-OPT/CROSS-OPT-1.0
   ```
3. To optimize a structure, execute the main file corresponding to your case:
   ```bash
   Main_TrussNbars_ObjFuncType_ConstraintType_Solver
   ```

This package includes the following files:
* Main_Truss10_MassMin_YieldStress_CE.m  -- Mass minimization with yield stress constraints for a 10 bars truss using CE
* Main_Truss10_MassMin_YieldStress_GA.m  -- Mass minimization with yield stress constraints for a 10 bars truss using GA
* Main_Truss10_MassMin_YieldStress_SPQ.m  -- Mass minimization with yield stress constraints for a 10 bars truss using SQP
* Main_Truss10_MassMin_FrequencyBounds_CE.m  -- Mass minimization with frequency constraints for a 10 bars truss using CE
* Main_Truss10_MassMin_FrequencyBounds_GA.m  -- Mass minimization with frequency constraints for a 10 bars truss using GA
* Main_Truss10_MassMin_FrequencyBounds_SQP.m  -- Mass minimization with frequency constraints for a 10 bars truss using SQP

### Documentation
The routines in **CROSS-OPT** are well-commented to explain their functionality. Each routine includes a description of its purpose, as well as inputs and outputs. 

### Reproducibility

Simulations done with **CROSS-OPT** are fully reproducible, as can be seen on this <a href="https://codeocean.com/capsule/xxx" target="_blank">CodeOcean capsule</a>.

### Authors
- Marcos Vinicius Issa
- Anderson Pereira
- Americo Cunha Jr

### Citing CROSS-OPT
We ask the code users to cite the following manuscript in any publications reporting work done with our code:
- *M. V. Issa, A. Pereira and A. Cunha Jr, The cross-entropy method for nonconvex structural optimization, 2024 (under review)*

```
@article{Issa2024CROSS-OPT,
   author  = {M. V. Issa and A. Pereira and A {Cunha~Jr}},
   title   = {The cross-entropy method for nonconvex structural optimization},
   journal = {Under Review},
   year    = {2024},
   volume  = {~},
   pages   = {~},
   doi    = {~},
}
```

### License
**CROSS-OPT** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%"> &nbsp; &nbsp; <img src="logo/logo_pucrio_color.jpg" width="07%">

### Funding

<img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">  &nbsp; &nbsp; &nbsp; <img src="logo/faperj.jpg" width="20%">
