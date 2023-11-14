# NEMO_activities

This repository includes subdirectories relating to my own personal projects in the NEMO code. The intention is for anyone to reproduce the work set out here. It is also still a work in progress.

The code in this repo relates to NEMO v4.0.4 from the [forge.ipsl directory](https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.4/). One can obtain the source files by using `svn co ...`. A more up-to-date NEMO version (v4.2) is available on [GitLab](https://forge.nemo-ocean.eu/nemo/nemo). At some point in the near future, the code will be updated for use in NEMO v4.2.

A brief outline is given below of what is included in each subdirectory. Further information is found in each subdirectory, with details given in each respective `documentation.md`.

**IDEAL**
- In `/IDEAL`, there are `src`, `EXP`, `tools`, and `scripts` for the idealised channel model configuration.
- The `/scripts` folder contains bash scripts for automating job submissions on Monsoon2, along with python scripts for setting up the model and also analysing model output.
- `/pyCDFTOOLS` contains specific functions to analyse model output. See [J. Mak GitHub repo](https://github.com/julianmak/NEMO-related/tree/master/pyCDFTOOLS) for additional `pyCDFTOOLS`.
- A `documentation.md` file also exists that documents the trials and errors in setting up this configuration.

n.b. a channel model confiuration `CANAL` already exists in the NEMO `cfgs` directory. We think this configuration (`IDEAL`) offers more flexibility and ease of use without someone having to dive into FORTRAN code to make changes.

**QG_Leith**
- Includes modifications to `ldfdyn.F90` for the implementation of the Leith schemes in NEMO. 
- Changes to `step.F90` for the inclusion of the Leith Schemes.
- Additions to `oce.F90` for daily stretching calculations in QG Leith.
- Modifications to `ldftra.F90` for the implementation of Leith schemes as GM/Redi coefficients. Currently in the testing stage.
- A `documentation.md` file that documents the progress.
The Leith schemes are still very much in the testing stage and may be prone to errors. Both schemes run well in idealised and realistic configurations without crashing.

**ORCA025**<br/>
The Leith schemes are being tested in realistic forced ocean sea-ice configurations.

Similarly to `/IDEAL`, we employ `pyCDFTOOLS` with minor modifications made. We hope to bring both versions together so they are cross compatible, and can be used with any configuration. Scripts are also being generated to analyse the data e.g. `calc_cdfmoy.py`.

A module called `nemo_toolkit`, with `data_view.py` to easily load and plot nemo data. Still under development but plotting capabilities already include map projections, transects, and time series. 

Modifications relating to QG Leith in orca025 are also made to:

- `ldftra.F90` which also includes a Southern Ocean package developed at the Met Office.

More details to be added here in due course.


