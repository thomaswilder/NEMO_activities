# NEMO_activities

This repository will include my activities in the NEMO ocean model.

**IDEAL**
- In `/IDEAL`, there are `src`, `EXP`, `tools`, and `scripts` for the idealised channel model configuration.
- The `/scripts` folder contains bash scripts for automating job submissions on Monsoon2, along with python scripts for setting up the model and also analysing model output.
- A `documentation.md` file also exists that documents the trials and errors in setting up this configuration.

**QG_Leith**
- Includes modifications to `ldfdyn.F90` for the implmentation of the Leith schemes in NEMO. 
- Changes to `step.F90` for the introduction of the Leith Schemes.
- A `documentation.md` file that documents my progress in coding this up.
The Leith schemes are still very much in the testing stage and are prone to errors!

**ORCA025**
