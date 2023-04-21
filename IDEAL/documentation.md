# Objective
Documenting the development of the idealised nemo configuration. The aim is to design a **Neverworld2** configuration, detailed in Marques et al. (2022) in order to test the implementation of QG Leith. We begin with a simple Southern Ocean re-entrant **channel model**.

### Tools
- Checkout NEMO 4.0.4 from either ipsl.forge or the MOSRS repository.
- We make use of the excellent [documentation](https://afstyles.github.io/nemo-compilation/) by Andrew Styles for setting up NEMO on Monsoon2.
- Julian Mak also provides a [guide](https://nemo-related.readthedocs.io/en/latest/nemo_notes/unagi_config.html) with great detail into setting up an idealised channel model.
- Gyre_Pisces example configuration in cfgs/ folder.
- tools/DOMAINcfg to build a domain_cfg.nc file.
- tools/REBUILD_NEMO to combine output files.
- tools/DMP_TOOLS to generate a restoring file.

### Neverworld specific
- Use a non-linear equation of state.

### (15/11/2022)
### Initial experimental design 
- A simple east-west periodic channel model to get things started.
- Domain size is (L_x, L_y) = (1000, 2000) km with grid spacing of 10 km, so (nx, ny) = (100, 200+2).
- In the vertical, L_z = 3000 m with nz = 14+1 z-levels. Vertical grid calculated in-house.
- Ocean bottom is flat (for now).
- Temperature, salinity, and wind stress are calculated using gen_nemo_ideal_fields.py, a modified J.Mak's _unagi_.py file.

### Modifying namelist_cfg file
- Biharmonic viscosity, `rn_Uv = 1` m/s, `rn_Lv = 10.e+3` gives `Am = 8*e+10 m^4/s`.
- Vertical diffusion `rn_avt0=0.` so temp does not diffuse anywhere.
- Time-step, `rn_dqt = 800` s, and `nn_itend=108` to run for one whole day initially.
- Using absolute wind stress to start, so `rn_vfac = 0` under `namsbc_blk`. This is set in the reference and only needs changing for relative wind stress.

### (17/11/2022)
### Parameters (namelist_cfg)
- Defined `ln_traldf_OFF   = .true.` in `namelist_cfg` as want an adiabatic system for now.
- Set `nn_itend    =   10800   !  last  time step. 10800=100 days at 800 s timestep`.
- Set `nn_write    =    1080   !  frequency of write in the output file`. 

### Making the model run
- Dropping J.Mak's `postprocess.sh` for now as a few things are not clear in the script. Using `rebuild_nemo file $NUM_CPU` in terminal instead to combine output data (see README.rst in REBUILD_NEMO). Will write a script later to make the combining of files smoother. Example file would be `IDEAL_1d_00010101_00010410_grid_U`
- There are 20 time dimensions in the output files, but 10800/1080 = 10...? So interval write is 5 days in output file `IDEAL....nc`.
- Looks like `nn_write` does not control when output is written to file... This is for the ocean.output file?

### (21/11/2022)
- Output frequency is managed in the file `file_def_nemo.xml`
- Set `enable=".FALSE." for all output except 1d... nope.
- Moved "<file id="file1"...` below `<file_group id="1d!..." hmm
- Added additional file (diagnostics) to 1 day output section.
- A few syntax errors, now ironed out, and model outputting 1 day mean data. See file `file_def_nemo.xml` in `EXPH10` folder. (insert link here).

### Spinning up the channel model
### (23/11/2022)
- Apply sea surface restoring and a northern sponge layer (see Munday and Zhai (2017) paper for details on values). In NEMO guide, see section 4.6 `Tracer Damping` for setting this up.

### (24/11/2022)
- Try surface restoring and see if that is enough. See section 6.2.1 (pg 70) titled `input data specification` for explanation of including restoring data files.
- We will restore yearly to begin with as that is what is coded in `gen_nemo_ideal_fields.py`. To choose monthly restoring that is repeated every year, create 12 time dimensions and have a yearly open/close of file.
- Run for 100 days and runs fine with restoring on.
- Set `nn_itend=38880` which is 360 days.
- Running for one year is very quick.

### (25/11/22)
- Written script to combine output files and move files about, `comout.sh`.
- Running channel model for 2 years. Needs to spun up for longer (~800 years), see Zhai and Munday (2014) for details.
- Need a northern boundary temperature sponge layer too. Calculating `mesh.nc` file using `ln_meshmask = .true.` in namelist `&namdom` in `namelist_cfg`.
- Started writing the northern boundary restoring sponge layer in nemo code. See file `north_bdy.F90` in `DMP_TOOLS` folder. ** need to add this ** . 

n.b. the decision has been made to get the channel model working properly and test the QG Leith code in this config as it is computationally inexpensive. Thereafter, we will build Neverworld2.

### (29/11/22)
- Calculating the sponge layer. Found a bug in `utils.F90` under function `CALL dimlen( ncin, 'z', jpk )`, which should use `nav_lev` for the `z` dimension. Changed this and now a resto.nc file is created. Need to confirm if the resto.nc data is correct. Need to also upload the modified `DMP_TOOLS/src` files.
- Run `./maketools -m XC_MONSOON_INTEL -n DMP_TOOLS` to generate a `make_dmp_file.exe`. This is then executed using the bash script `mkmeshscript.sh` that produces a `resto.nc` file. See folder `shell_scripts` for this bash script.

### (30/11/22)
- Setting the relaxation at the northern boundary, have included `namtra_dmp` namelist with `ln_tradmp = .true.` and `resto.nc`. Also have to set `ln_tsd_dmp = .true.` in `namtsd` namelist.
- Want netcdf to output data yearly/twice yearly during the spin up stage.
- Test if the model runs first... which it does, and also from a `restart.nc` file.
- Spin up the model for 800 years with yrly output, then restart and run for 1 yr with 5 day output to examine the flow and tracer field...
- Looks like this needs to be done in batches, so run 100 yrs then another 100 yrs and so on. 
- Having issues with the restart counter...

### (5/12/22)
- Model is running for 10 years.
- Writing shell script to set up a PBS job chaining for the spin up stage. We want to continuously submit nemo jobs to achieve a 800 year spin up stage without having to do this manually. See [this](https://waterprogramming.wordpress.com/2016/08/28/pbs-job-chaining/) webpagefor details on pbs chaining.

### (6/12/22)
- If `nn_stock=0` then `rebuild_nemo` does not work. Needs to be set to the number of iterations for the model to run.

### (7/12/22)
The shell scripts in [this folder](shell_scripts) will allow the running of multiple nemo jobs in Monsoon2. The order of jobs can be found in `submit_N_nemo.sh`. The idea behind this is to enable an 800 year spin up within the constraints of a 4 hour wallclock limit on Monsoon2. Shorter running jobs are set with a higher priority.
- Have set up 8 jobs on Monsoon2 each running 50 years in length for a total of 400 years...

### (8/12/22)
- Model didn't seem to start properly. Lets try again with smaller runs and check output at each restart e.g. 3240 time steps. This works.
- Might need to run a shorter simulation and more jobs. Try 10 months with 10 pickups. Look at the output.
- Running 20 simulations of 20 years. **Holding breath**

### (12/12/22)
- 420 years of spin up has been completed, and the profiles look good. Run for a further 380 years.

### (14/12/22)
The 800 year model spin up is complete. There looks to be some vertical velocity noise near the northern and southern boundaries, this could explain why J. Mak added additional damping at northern boundary. The overall kinetic energy pattern looks good however. Transects of temperature look okay too, do need to compare initial with 800 year point however. As this is only a test case, it doesn't seem important yet to tackle this slight numerical noise.

### Adding in the QG Leith code
- We add the new `ldfdyn.F90` to the `src/OCE/LDF` directory in the NEMO v4.0.4 model. We rename the original `ldfdyn.F90.orig` for backup purposes.

- Looks like relative vorticity `relvor` is missing from `diawri.F90` in the `src/OCE` directory, but is present in the `CANAL` test case `diawri.F90`, so have put the `CANAL` version into the `MY_SRC` directory in `IDEAL. Have also recompiled nemo using `./makenemo` but this doesn't seem to be working. Can calculate relative vorticity offline but seems time consuming when it is calculated online.
- Outputting `ahmt_3d` is only possible for viscosity cases 31 and 32 in `ldfdyn.F90`... which makes sense.
- Maybe add `diawri.F90` into the main `src` directory? Recompile this.

### 21/12/22
- Haven't yet tested QG Leith, trying to output correct diagnostics. Try running model with a 3D viscous scheme to check `ahm` outputs.
- There were some errors in the `diawri.F90` code. Spaces between code lines should be separated by a `!`, so adding these in solved the issue.
- All diagnostics working. Relative vorticity field looks great too.

Made various edits to `ldfdyn.F90` and `step.F90` to incorporate QG Leith. Errors followed e.g. `*** glibc detected *** ./nemo: free(): invalid pointer`

### 9/1/23
- Changed the allocation of `zbu` ... and moved them to the start of the module where ddivmagsq is allocated. There was also an index error an allocation where `ji` should have been `jpi`. Lets see if these do the trick... partially.
- The model runs, accepting the modified `ldfdyn.F90` with QG Leith, however the output of viscosity is not physical. There are grid lines and high levels of viscosity at the boundaries, with the majority equal to zero. There are some instances of eddy activity in the pattern though.

### 10/1/23
- Checking if a 2D Leith viscosity scheme will run. 2D Leith runs and gives a physical output. 
- Need to make a 2D Leith scheme for user selection!

### 11/1/23
- `ldfdyn.F90` updated to include 2D and QG leith as options to the user. Have made some changes to the naming conventions of terms, e.g. `zvsdx` is now `zwzsdx`, since vorticity is `zwz`.
- Possible error detected in QG routine for the Rossby number calculation, `zrosq`. Was including `e1e2f**2` which is an area squared, and only need area. Now the code blows up... so what is the best way to debug the code and find where the error is? 
- In `ocean.output_0004`, it shows that the model blows up after 400 timesteps, and near the surface. 

### 17/1/23
- 2D Leith scheme working when user chooses `case 33` and `rn_c2dc = 1`. Think there is some model noise and this could be due to a lack of model diffusion on tracers?
- In QG Leith routine, modified lower layer boundary condition for stretching term calculation, and also replace `zst` with `zstlim` around line 698 due to overwriting data which should not have been done. Could this be the cause of the grid lines? Nope.
- QG Leith now runs but there is large model noise and viscosities (order magnitude bigger than it should be), particularly at the boundary walls. There are also some 'grid lines' in the visc output at the northern boundary throughout water column, in the `ji` and `jj` direction. Could this be due to the restoring boundary conditions in `IDEAL`? Or, the grid lines could be where the QG dynamics are failing and the scaling is coming in to effect, so an error must be present in that part of the code?
- Looks like the Froude number is quite big everywhere and the scaled stretching is being used where it shouldn't be. But the QG Stretching term is also too big?
- In `ocean.output`, buoyancy is periodically zero when `jj=1`, and this could be the Southern boundary wall?
- Froude number does reach values below 1, see line 1652748 in `ocean.output`, but Rossby number is still smaller so stratification dominates.

### 18/1/23
In an attempt to fix the Halo issue in the QG Leith field, we:
- Added in `lbc_lnk_multi` function for `zstlim` on F point after calculation of Rossby number, but it seems like the Halo has grown horizontally.
- Bringing this `lbc_lnk_multi` function forward in the routine for `zst` has removed the vertical Halo lines, but one still persists at the northern boundary.

### 19/1/23
### new spin up run with diffusion
- Slight wave like patterns/ grid scale noise which may be a result of turning off diffusion in the model. Run another 800 year spin up with diffusion.
- Choosing parameters: $\kappa_z = 5.5 \times 10^{-5} \text{ m}^2 \text{ s}^{-1}$ and $\kappa_4 \approx 3 \times 10^9 \text{ m}^4 \text{ s}^{-1}$.
- In `namtra_ldf`, `ln_traldf_blp = .true.`, `nn_aht_ijk_t = 30`, `rn_Ud = 0.036`, and `rn_Ld = 10.e+3`. These choices give us the above lateral $\kappa_4$. The scaling values are not indicative of the flow in the model.
- Run for 200 years to start then check the output before completing the 800 years.

### 6/2/23
- Diffusion run looks good, at least with a monthly output. See '.gif' in local 'IDEAL/Figures'.
- Run a daily mean simulation beyond 800 years for an additional month and check for any grid scale noise. There definitely looks to be less visible noise, perhaps none at all.

Let's attempt to spin up the channel model using the 2D Leith scheme, and include the diffusion term. Running for 200 years to start.

### 16/2/23
- 2D Leith spin up seems to be very noisy at day 200.
- Have run the 800 year spin up with diffusion... looks good at the monthly mean output.

The QG Leith routine is not outputting reasonable Froude number values (too large). The QG stretching term is also too big. A possible correction to this is to run a simulation with more vertical levels to improve the representation of the vertical baroclinic modes.
- Setting up a configuration called 'IDEAL/EXPZ45' that has 45 vertical z levels. Creating new domcfg.nc, forcing, state, resto files.
- Check 'testing.output' for model runtime length. Also, increase number of cores/nodes used to 108 for nemo.

2D Leith after 800 year spin up in 'IDEAL/EXPH10' seems noisy near surface and towards southern domain. Could this be due to model vertical resolution as well? Could the noise near southern boundary be due to the poor representation of eddies at the 10 km grid scale?

- Running EXPZ45 for a 200 year spin up. Looks good.

### 20/2/23
- 200 year spin up so far looks good. Lots of noise near northern boundary though.
- Running an additional 200 years because Monsoon has a job submission limit.

### Adding in diagnostics
- See [this thread](https://nemo-ocean.discourse.group/t/creating-nemo-output-diagnostics/446) on the nemo discourse forum about implementing diagnostics. Struggling to get my head around where terms are called, and just the overall method for achieving this.

### 22/02/23
- Density diagnostic is working. It is called in `step.F90` when `ldf_dyn` is called.
- Now attempt to output Froude and Rossby numbers.

### 28/2/23
- IDEAL with 45 z-levels is crashing when either Leith scheme is chosen. This is restarted from 800 year spin up.
- There are extremely high values of hdivn gradients in the upper left corner of domain, but then zero elsewhere. 
- Looks to be running now after using `CALL dyn_hor( kt )` in 2D Leith routine.

### Static instability?

### 9/3/23
- IDEAL is displaying negative $N^2$ near the southern boundary, so is statically unstable.
- It was suggested by Till to change the southern boundary temperature from near 0c to 4c. The papers that use this setup (e.g. Abernathey 2011 and Munday 2017) do not do this. They have KPP and surface fluxes, so I wonder if this is what's missing here...?
- Now, the convective adjustment looks to be working on tracers and momentum, `nn_evdm = 1`. Change this to work only on tracers, so `nn_evdm = 0`. Run this new spin up in `cfgs/IDEAL/EXPH10`. Also, J. Mak increased `rn_dqdt = -80` to have a stronger surface restoring. We will do this.

### 10/3/23
- Doing the above changes still results in negative N^2.

New temperature profile `IDEAL/statev2.nc` with 2c at southern boundary and linear increase to 16c at northern boundary. Removing noise too as this may have an effect?
- Running IDEAL for 50 years with new temp field. This has not worked.

Lets try a different convective adjustment scheme.
- Trying non-penetrative convection scheme `ln_zdfnpc=.true.` in `&namzdf`, and commenting out enhanced vertical diffusion. We use the default values for npc scheme. We also put `nn_etau=1` instead of `0` in `&namzdf_tke` so that TKE penetrates below the mixed layer.
- Running for 50 years ...
- Looking at year 20, at least compared to the last run, it looks like $N^2$ is not negative.

Spinning up for 400 years. Exciting.

### 14/3/23
### Channel Model Bump

Extending IDEAL to include a bump in the bathymetry. Setup now has $L_x=6000$ km, bump is 1500 m in height and 2000 km in length. This is found in `cfgs/IDEAL/EXPH10`
- Model takes between 2.2 - 3 minutes to run one month (~30 days) using biharmonic viscosity, so is estimated at 88 - 120 hrs for 200 years of spin up.
- Spinning up the grid-dependent biharmonic viscous simulation with 200 years to start.

### 15/3/23
- Model output empty in all files.
- Running monthly simulations for a few months provides output.

Run for a year and the output is fine. Running for 50 years now.

### 16/3/23
- Year 43 looks fine so far. Eddy field is developing nicely. Year 50 is good.
- Running a further 50 years at 5 year pick up intervals.
- Renamed variables in `file_def_nemo.xml` from 50 years because CDFTools uses set variable naming conventions. Tried renaming `uvel` to `vozocrtx` in xarray but something messed up. Lets try the cdftools `cdfmoy` and `cdfeke` from year 100 of spin up...

### 17/3/23
- The model starts outputting NaNs in velocity field around year 55-60. The ssh field does look quite large in places. The model continues to run though, and temperature field is still present and evolving. 
- Wonder if something has gone wrong in xios? Yep, looks like `rebuild_nemo.exe` was terminated after using 983356 kB of memeory. So .nc files are too large? Allow the model to complete this stage then restart with yearly pickups again.

Working fine now at yearly output. Keep running until year 400 at least.

Running for an 70 or so years... maxed out job submission on Monsoon2! 150 years of spin looks good.

### 20/3/23
- Spin up at year 217 looks good.
- Increasing the number of cores from 180 to 240, see how this changes the spin up time. Might increase to 360 too. 

### 21/3/23
- Spin up at year 280 looking good, a nice big eddy has formed off the eastern bump.
- Issues in restart scripts? Error in `postproc.error` so the option `depend=afterok` would not allow the job to submit. Issue resolved by removing `rm submit*.pbs` at end of `jobnamup.pbs`.

### 23/3/23
- 400 years looks good, compares well with [Munday et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014PA002675) so far. 
- Reduced wall time limit to 30 mins per job.

### 27/3/23
- 600 years of spin up using biharmonic viscosity complete.

Checking Leith schemes with 600 year pick up.
- Will run each scheme for a further 20 years, including biharmonic, then compute diagnostics such as EKE.
- ISSUES with QG Leith SCHEME! Really large values in places making model blow up.

### 28/3/23
- NEMO not compiling. Errors such as
  ```
    ftn-2105 crayftn: ERROR in command line
    "-i" is an invalid command-line option.
    ftn-2191 crayftn: ERROR in command line
    "-size" is an invalid argument to the "-r" option.
  ```
- Have emailed Monsoon support...


### 20/4/23
Model working now...

- There is some noise in the QG Leith run. We will try a 4th order tracer advection scheme by setting `nn_cen_h = 4` in `namtra_adv` in `namelist_cfg`.
- Have also reduced number of processes needed by model from 360 to 180 to lower computational cost.

Could do with an exit with error code line in bash script `combine_tidy.sh`. Running this function without correct XIOS input deletes diagnostic outputs regardless if rebuild was successful or not. 

### 21/4/23
- 4th order tracer advection scheme doesn't reduce the noise.
- It does look like the noise is present at the base of the mixed layer.
- Try a smaller timestep e.g. 200 s. Also didn't solve the problem of noise.
