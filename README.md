# NEMO_activities

This repository includes subdirectories relating to my own personal projects in the NEMO code. The intention is for anyone to reproduce the work set out here. It is also still a work in progress.

The code in this repo refers to NEMO v4.0.4 from the [forge.ipsl directory](https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.4/). One can obtain the source files by using `svn co ...`. A more up-to-date NEMO version (v4.2) is available on [GitLab](https://forge.nemo-ocean.eu/nemo/nemo).

A brief outline is given below of what is included in each subdirectory. Further information is found in each subdirectory, withd details given in each respective `documentation.md`.

**IDEAL**
- In `/IDEAL`, there are `src`, `EXP`, `tools`, and `scripts` for the idealised channel model configuration.
- The `/scripts` folder contains bash scripts for automating job submissions on Monsoon2, along with python scripts for setting up the model and also analysing model output.
- `/pyCDFTOOLS` contains specific functions to analyse model output. See [J. Mak GitHub repo](https://github.com/julianmak/NEMO-related/tree/master/pyCDFTOOLS) for additional `pyCDFTOOLS`.
- A `documentation.md` file also exists that documents the trials and errors in setting up this configuration.

**QG_Leith**
- Includes modifications to `ldfdyn.F90` for the implmentation of the Leith schemes in NEMO. 
- Changes to `step.F90` for the introduction of the Leith Schemes.
- A `documentation.md` file that documents my progress in coding this up.
The Leith schemes are still very much in the testing stage and are prone to errors!

**ORCA025**
Work will begin soon...
