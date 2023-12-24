# NEMO_activities

## Description

This repository includes work relating to the development of the NEMO model and any code for analysing NEMO data. Seting out the code here provides total transparency of the NEMO model development. The intention is to also provide some tools for the analysis of NEMO data. 

The code in this repository relates to NEMO v4.0.4 from the [forge.ipsl directory](https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.4/). One can obtain the source files by using `svn co ...`. NEMO version 4.2 is available on [GitLab](https://forge.nemo-ocean.eu/nemo/nemo). At some point in the near future, the code will be updated for use in NEMO v4.2.

A brief outline is given below of what can be found in each subdirectory. 

## Directories

**IDEAL**
An idealised channel model representing the Southern Ocean.

- In [IDEAL](IDEAL/), there are `src`, `EXP`, `tools`, and `scripts` for the idealised channel model configuration.
- The [scripts](IDEAL/scripts) folder contains bash scripts for automating job submissions on Monsoon2, along with python scripts for setting up the model and also analysing model output.
- [pyCDFTOOLS](IDEAL/pyCDFTOOLS) contains specific functions to analyse model output. See [J. Mak GitHub repo](https://github.com/julianmak/NEMO-related/tree/master/pyCDFTOOLS) for additional `pyCDFTOOLS`. n.b. differs slightly to the version in ORCA025.
- A `documentation.md` file also exists that documents the trials and errors in setting up the channel model.

n.b. a channel model confiuration `CANAL` already exists in the NEMO `cfgs` directory. We think this configuration (`IDEAL`) offers more flexibility and ease of use without having to dive into FORTRAN code to make changes.

**QG_Leith**

Code changes relating to the implemtation of Quasi-Geostrophic Leith viscosity in NEMO.

- Includes modifications to `ldfdyn.F90` for the implementation of the Leith schemes in NEMO. 
- Changes to `step.F90` for the inclusion of the Leith Schemes.
- Additions to `oce.F90` for daily stretching calculations in QG Leith.
- Modifications to `ldftra.F90` for the implementation of Leith schemes as GM/Redi coefficients. Currently in the testing stage.
- A `documentation.md` file that documents the progress.

The Leith schemes have been tested and are currently performing well.

**ORCA025**<br/>

The Leith schemes are being tested in realistic forced ocean sea-ice configurations.

Similarly to `/IDEAL`, we employ [pyCDFTOOLS](ORCA025/pyCDFTOOLS) with minor modifications made to the idealised version. We hope to bring both versions together so they are cross compatible. 

Python scripts are available in [python_related](ORCA025/python_related) that calculate and plot NEMO data. These scripts use pyCDFTOOLS and [nemo_toolkit](ORCA025/nemo_toolkit), where both are still in active development.

Modifications relating to QG Leith in orca025 are also made to:

- `ldftra.F90` which also includes a Southern Ocean package developed at the Met Office.



