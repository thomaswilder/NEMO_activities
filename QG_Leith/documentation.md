## Documenting QG Leith Implementation

This documents the development: 
1) QG leith viscosity parameterisation,
2) Using Leith coefficients in the Gent-McWilliams and Redi scheme,

The QG Leith viscosity implementation makes use of MITgcm implementation and work by Bachman et al. (2017) and Pearson et al. (2017).

See also the `documentation.md` files in `IDEAL` folder [here](https://github.com/thomaswilder/nemo-IDEAL) and `ORCA025` folder for further details on model setup to test the Leith schemes.

### QG Leith viscosity
### 20/02/23
- Need to figure out how to output diagnostics for terms in the QG Leith routine to establish why it is not working. Spatial maps will help greatly. Have asked this question to the nemo discourse community.
- The routine works when modifying `step.F90` and calling `bn2` and `eos` so `ldf_dyn` inputs density and buoyancy frequency into the QG Leith routine. However, some other routines e.g. `ldftra.F90` include `rn2` without inputting or calling buoyancy frequency squared from `step.F90`. So need to figure this out.

### 21/02/23
- To output buoyancy, Froude and Rossby number diagnostics, need to add an IF statement at the end of `ldfdyn.F90` so only valid when `CASE( 34 )` is chosen. Unsure yet how to achieve this.
- Check if density outputs...
- Also checking if divergence is working too, outputting individual values in routine into `ocean.output`

Why are the Froude and Rossby numbers still being output in ascii ocean.output file when they are not in the `ldfdyn.F90` file. Hmmmmm?

### 22/02/23
- Did not use the correct `nemo` executable!
- `hdivn` is populated since there are gradients of divergence coming through.
- Density is now being output in the diagnostic file since it is called in `step.F90`.

- huge E R R O R and model is blowing up with QG Leith in current configuration. Will be interesting to test it on the higher vertical resolution config.
- Blows up with both before and now density and buoyancy frequency.

### 27/02/23
NEMO IDEAL_EXPZ45 is blowing up. Test on 15 vertical level run in IDEAL2.

### 28/2/23
- Issues with `hdivn` in Leith routines. Probably due to the order the routines are called in `step.F90`. 
- Solved this by adding `CALL dyn_hor( kt )` in the Leith routines, which now populates `hdivn`.

Wondering whether the issue resides in taking the minimum of the stretching term and dimensionless numbers before taking the horizontal gradients. Nevertheless, the QG scheme should default to a 2D Leith scheme.

- Currently the QG scheme does not run and IDEAL exits with code 139 and segmentation fault. Errors in accessing memory...?

### 2/3/23
- Fixed previous issue as `rro2` and `rfr2` were not allocated memory space.

- Outputting diagnostics every timestep since model is crashing after 12 or so timesteps.

**Implementing Pearson et al configuration!**

- See `ldfdyn.F90.Pearson`.
- Including mixed layer depth index is found through `CALL zdf_phy( kstp )` in `step.F90` which calls `CALL zdf_mxl( kt )` within, generating mixed layer depth index `nmln`. Not done this yet.

- Running without mixed layer depth condition: $Bu$ appears as negative in some regions due to negative $N^2$, which implies there is static instability. This is in 45 vertical level model.
- Checking if this happens in 15 vertical level setup. Yes.

Some issues:

1. Some indices are not outputting in the `ocean.output` file. Stretching outputs at `(ji,jj,jk)=(1, 14, 1)`, but terms such as `Ro^2` and `N^2` are not given at this point. Is this a `lbc_lnk_multi` function issue?
2. There are regions of static instability (heavy waters over light), this is the case in both simulations. Primarily at the southern boundary. Some at the northern boundary in EXPZ45 near the surface.

### 3/3/23
- The diagnostic field for Rossby and Burger numbers look well populated, and the values look reasonable too. 
- The QG stretching field appears to be the one with issues. There are small blocks of very large values.
- But I still don't understand why QG does not transition to 2D Leith?
- Wonder if the issue reside in not filling the buoyancy field with the boundary conditions. Lets try that. Seems to have fixed the blocks of very large stretching.
- QG Leith displays very large values when `N^2` is negative, and these negative values seem to be impacting the stretching term at the southern boundary. What is the fix for this?

*Besides, lets implement the mixed layer condition that Pearson et al. did.*

- Working okay. Need to include meridional gradient $\beta$ in vorticity. 
- Large values of viscosity near southern boundary. 

### 6/3/23
- Issues seem to be located in the stretching term... lets look at the buoyancy gradients to isolate where the issue is.
- Buoyancy gradients appear to contribute to the profile of the QG stretching term. Looks like a model issue with IDEAL setup...

### 9/3/23
- Going to add a condition into QG Leith routine that reverts to 2D Leith when $N^2<0$.
- But we need a buoyancy field everywhere for the gradients? So calculate $b$ in the domain. This has solved some of the large gradients near the southern boundary. Large gradients still persist.

### 13/3/23
- Having solved the negative $N^2$ in `IDEAL`, no longer need to random $N^2<0$ condition as suggested before. In any case, this caused incorrect values in domain for stretching terms.

The QG Leith looks okay for a few timesteps.

Output `zstlimx` and `zstlimy` to compare with `zstx` and `zsty`, asking whether the scheme is using the scaling or the actual stretching term.

Extremely high values of QG stretching just below mixed layer, which feeds in to QG viscous coefficient. There are also high values near southern boundary at depth, which dominate over vorticity gradients. This might be reasonable considering Pearson et al. (2017), though spatially it looks slightly wrong.

Wonder if should calculate buoyancy gradients everywhere too? Checking ...


### 28/3/23
- Experiencing issuses with QG Leith again. Seeing extremely large values in the stretching term even though the algorithm should be choosing the smallest value. The pattern in stretching comes from $N^2$.
- Wonder if there is an issue on line 815 and 820 with the use of min? Have changed to MIN to investigate.
- Plotting in matlab `figure();contourf(zstlimx(290:330,5:40,3,1)');colorbar;title('zstlimx')` we can examine the fields.

### 30/3/23
Compile issues from a clean start:
- NEMO won't compile from scratch using modified `step.F90` and `ldfdyn.F90`, but will compile from original files.
- Updating these files in `src/OCE` then compiling without cleaning results in a successful build.

```
/home/d02/twilder/NEMO/NEMO_4.0.4_mirror/cfgs/IDEAL/BLD/ppsrc/nemo/trdvor.f90(22): error #7002: Error in opening the compiled module file.  Check INCLUDE paths.   [LDFDYN]
   USE ldfdyn          ! ocean active tracers: lateral physics
-------^
```

- Wondering if the PUBLIC variables are the issue? No.

There are no errors in the code otherwise that would get flagged up?

Start building `ldfdyn.F90` from scratch...? See `NEMO/QGLeith/debug/` directory. Also start again in `step.F90`.

Compiling with 2D Leith:
- Fails with 
```
/home/d02/twilder/NEMO/NEMO_4.0.4_mirror/cfgs/IDEAL/BLD/lib/lib__fcm__nemo.a(ldfdyn.o): In function `ldfdyn_mp_ldf_dyn_':
ldfdyn.f90:(.text+0x45fc): undefined reference to `div_hor_'
fcm_internal load failed (256)
gmake: *** [nemo.exe] Error 1
```
- Remove `CALL div_hor( kt )`. No. Added in `USE divhor` and put call back in subroutine.
- This resulted in the error relating to `trdvor.f90`. Why?

Trying a fresh compile using `ldfdyn.F90.Pearsonmxl` without the call for divergence. **SUCCESS!**

Now, does adding in `MIN` rather than `min` make a difference?

Restarting file not working? Trying with biharmonic viscosity. This Works! 

Fails with each Leith scheme. Is this because of zero horizontal divergence? Could i call this in `step.F90` and feed in to the ldfdyn subroutine.

### 31st March

Adding in calculation of `hdivn` to `ldfdyn.F90` instead of trying to call the function. Making it PRIVATE. No, changed name to `hdivnqg`.

IDEAL now completed one timestep with 2D Leith. Model numerically unstable without `hdivn`.
IDEAL now completed one timestep with QG Leith.

Examining the QG Leith output in MATLAB:
- Some of the values are chosen correctly when assessing the QG limit around line 837 in `ldfdyn.F90`. E.g. `zstx` is smaller than `tmpzstx = zwzdx/MAX(rbu,rro2)` and so matches `zstlimx`.
- But some of the values are incorrect, with neither `zstx` or `tmpzstx` matching the output of `zstlimx`, and so contributes to the large viscosity values.
- Diagnose `tmpzstx` in the routine and output.

Model is not running...?

### 3rd April
Model not running because used incorrect indices in `tmpzstx` allocation in `ldfdyn.F90`

The code is working, but there are patches where the stretching values are just too large.
- Adding in Coriolis parameter to QG vorticity.

What is the mixed layer depth? Output the index. Use `hmlpt` instead, which is the depth of last T-point inside the mixed layer [m].
- Adding `mld_dt02` to `file_def_nemo.xml` to output mixed layer depth.

Using `mld_dt02` to define the mixed layer depth in `ldfdyn.F90`. No, the use of `nmln` is correct.

The scheme that Pearson implemented does not smoothly transition from the QG to 2D regime. Recall we had issues with the Bachman scheme.

So:
- Really large vertical vorticity gradients and stretching at the 2D/QG transition boundary. Difficult to understand why?


### 4th April
Will remove `ff_f` and see if $partial_x q_2$ is the same? And also IF ELSE Statement around line 781, to check if strange lines are still present.
- Lines gone.
Vorticity should not have a condition in, so  $q_2 = f + \zeta$ for all space.

On examination, looks like the Burger number needs a minimum $N^2$? Do this.

But also, looks like QG scheme is calculating stretching when it shouldn't be. E.g. at depth 1000 m (z=20) grid points 310:322, 15:22, when it should be in the 2D Leith regime. We are outputting `mld_dt02` which is the mixed layer depending on temperature. Probably don't want to rely on this form. Instead:
- Output different mixed layer depth diagnostics e.g. `mldr10_3`, `mldr0_1`.
- Possibly try these as a new condition in the QG Leith scheme, depending on above.
- See `src/OCE/DIA/diahth.F90` for MLD calculation.


### 5th April
Trying to figure out how to output the vertical level of the mixed layer depth.
- Editing `zdfmxl.F90` to output `nmln`? Can an integer be used in `iom_put`?

Can we assign the depth at level of mixed layer and output that as a diagnostic?

Have computed QG mixed layer depth in `ldfdyn.F90` and reverted `zdfmxl.F90` back to its original. Check model output... 
- `mldr10_3` is deeper than `mldt02` in places.
- Where the QG Leith scheme fails is due to the depth of the mixed layer, which is being chosen incorrectly.

Outputting `mldr10_1` to see if it compares with `mld_qg`. Yes it does. 

So, going to modify `zdfmxl.F90` to output a mixed layer depth with `rho_c = 0.03`. A slightly larger mixed layer criterion. 
- Added in `nmln3` and made it public, but compile error saying `nmln3` must have an explicit type...

Posted a thread on NEMO discourse about this issue.

### 12th April
Lets check if the new mixed layer depth value actually works by modifying `rho_c=0.03`. Uploaded alternative `zdfmxl` and update `ldfdyn` accordingly.
- Recompiling NEMO... good.
- Now, running for 2 timesteps and examining QG Leith output... Mixed layer depth still not deep enough. 

Is there anyway to incorporate `mldr10_3` into `ldfdyn`? This produces a deeper mixed layer.
- Modifying `rho_c` in `zdfmxl` then recompiling doesn't carry through to the `/Work` directory...?

Write in `mldr10_3` code into QG Leith and use as mixed layer depth condition. Go through `diahth.F90` and add into `ldfdyn`.

Written in a mixed layer depth calculation into the QG Leith routine... fingers crossed. 

Model blows up around around depth of bump.
- In the IF statements for below mixed layer, have changed the condition `jk < jpkm1` to `jk < ( mbkt(ji,jj) - 1 )`, since `mbkt` gives the index of bottom last T point in ocean, compared with `jpkm1`, which is the absolute bottom `k` index.
- Successful run... 
	- Improvement in some regions index range (467:475,21:25,5,1). 

### 13th April
Routine is outputting `zstlimx` even though its in the mixed layer. So, what values are `nmlnqg`?
- Changing `mld_qg` output to `nmlnqg`.

The `nmln` computation in `zdfmxl.F90` is not sufficient. It initializes `nmln` to the w-point below 10 m, but then only changes its value when the mixed layer criteria is met. But what if it is never met? Then the initial value is incorrect. It should be initialized with the number of vertical index points, `mbkt`. 
- Lets try this...

The mixed layer depth and condition is now correct, but still issues persist.

I think the issue lies in the value of the Burger/Rossby numbers. A small Burger number implies the system is dominated by rotation, not stratification. The Rossby number does reach values of 1.

Boundary layer mixing scheme? No.

Looking at the smooth transition in Bachman that uses the Froude number.
- Examining the values in MATLAB suggests the Froude number could improve the scaled stretching term. 
- Running simulation with Froude number ... much improved, with maximum viscous coefficient of around 1500 m^2/s. Still seems big. 

Lets try and run the model for 1 day. 

### 14th April
But I set the diagnostics to 5 day outputs... so re-running the simulation for 10 days. The model did not output any errors previously.
- Data has been written to file... promising. YES. VALUES LOOK GREAT!
	- Max value at first 5 day mean is around 700 m^2/s, and next time step is 1000 m^2/s. Model not blowing up.

Lets try and run for a year with monthly mean output.
- Runs fine without errors.
- Viscosity output is noisy.

### 17th April
Checked the 2D Leith output at year 20, and this looks fine. Running this for further 80 years.

Why is QG Leith noisy?
- Outputting gradients of vorticity and divergence to see which terms contribute to the 'noise'.

### 18th April
In 2D Leith spin up at year 50, viscous coefficient shows no sign of noise.

Try spinning up with QG Leith. Maybe there is something strange going on with the change in viscous regime?
- Running QG Leith for 20 years.

Does horizontal divergence need `lbc_lnk_multi`? Lets try it.

**recompile IDEAL and IDEAL2** 

No change. Noise present in QG Leith.

### 19th April
Both Leith spin up runs have completed 20 years.
- 2D Leith has no noise, even after 20 years.
- QG Leith has **noise**.

Going to output contributions of QG Leith e.g. after grid averaging on to T-point.
- `ahmt_qg` and `ahmt_div`.
Pickup the QG Leith run at year 20 and run for a further year.

The noise is present in `ahmt_qg`, not the divergence part, which is opposite to what we thought.
- Computing `ahmt_qg` offline using `zwzd..` does not give the noisy pattern. But also, `mld_qg` was not showing integer values.
- Have changed the definition of `nmlnqg` from `REAL(wp)` to `INTEGER` in `ldfdyn.F90`.

Running the year again but outputting `zstlimx` and `zstlimy` as well to see if they are populated and impacting `ahmt_qg`.

Something odd is going on with `mld_qg`. Lets check `nmlnqg`...
- `nmlnqg` was `REAL(wp)` because `iom_put` won't output `INTEGER`.
Fixed this. Just assigned it prior to it being computed in routine.

Some weird numbers going on. Computing `zstlimx` in matlab is not consistent with output from model.

But, may be able to remove some `lbc_lnk` functions by using `jpi` instead of `jpim1 = jpi - 1`... Lets try this and see what comes of it.
- For example, changed this in `zbu`, `zstx`, `zwzdx`, Burger number routine, stretching ...
This compiled successfully!

Check model run in am.

### 20th April
Wonder if the monthly means are impacting the matlab computations? Lets restart and run for 2 timesetps.
- Yep, the monthly averaging is the reason. Single timesteps work.

Could try negative 1 in the `lbc_lnk` function for buoyancy gradients on U and V points. Maybe I am telling it to communicate in the wrong direction?

Going to try running the model with 4th order advection scheme for tracers.

### 6th June
Things possibly to consider:
- Smooth mixed layer depth somehow e.g. temporal or spatial average.
- Add a max viscosity value that could be user defined and informed by horizontal resolution.

### 27th June
Have not considered how to deal with horizontal gradients at/near bathymetry. This may be causing unrealistic values in the viscous coefficient.

Need to make use of module `zpshde.F90` which deals with the interpolation and gradients of tracers for z-coordinate with partial steps.???

Also, `IF( jk > nmlnqg(ji,jj) .AND. jk < ( mbkt(ji,jj) - 1 ) ) THEN` should be `jk < mbkt(ji,jj)` since there are two points at the bottom computed using 2D Leith rather than one.

### 28th June
Making modifications to `ldfdyn.F90`:
- Adding in masks to computations, which include `zbu` and its horizontal gradients, and also `zwz` gradients. See Git history for changes.
- QG Leith is computed for levels to `jk<mbkt`, not `mbkt-1`.
- Changed scale factors in vorticity gradient from `f` to `u` and `v`. Simiarly for the gradients of divergence and buoyancy. Incorrect scale factor in vertical stretching term, replaced `e3w_n` with `e3t_n`.

### 29th June
In order to compute stretching term daily, we create a subroutine for stretching in `ldfdyn.F90`, and compute it in the `ldfdyn` routine if the counter has stepped through one day.
- How do I save the stretching term? Check out diagnostic routines.

- Create a variable for stretching in `oce.F90` e.g. `zstlimx_n`
- Need to consider how a restart will work, possibly use `imo_rstput`? Or will this just recalculate on every restart and then use daily intervals from the restart date? 
- Use an IF function to compute stretching at right timestep. Also, use IF function in `ldfdyn_init` when allocating arrays. Use `kit000` for first timestep, then check for even timesteps of current timestep over timesteps per day e.g. `kt/108` and plus one so beginning of day, using `IF (MOD(kt,108)=0) THEN` if equal to zero then even and compute stretching term again and store it.

See progress so far in new `ldfdyn.F90`... making new subroutine called `ldf_dyn_str ()`

### 30th June
Progress so far:
- Model looked like it blew up on first time step. Couldn't figure out why...
- Instead, let stretching terms be zero for first day, then compute stretching terms daily thereafter. Model run for 3 days, and stretching terms stay the same throughout the day.

Lets run for one year. Blew up after about 16000 timesteps!
- Running for 4 months (around 12960 timesteps)

But there needs to be a condition that asks if the model is from a restart.
- Can we write stretching term to restart...?

But unsure if daily stretching values are appropriate due to model blowing up. What do the monthly outputs look like?

### 3rd July
Trying a shorter time for stretching calculations -> half day. Fails too on rose suite `u-xc856`.

It could be that the scaling term is computed at every timestep with evolving vorticity gradients, but only the QG stretching term including N^2^ is computed daily.
- Try an alternative routine that tries this.

New `ldfdyn.F90` routine not writing `zstx` to output file. Perhaps need to initialise the array with zeros, so doing this before the first call of `ldf_dyn_str ()`
- This fixed it!

`IDEAL_bump` runs for one year with output.
- Will commit this to QG Leith git repo.


### 6th July
I wonder if the use of `prd` and `bn2` at the before timestep are causing the numerical issues in ORCA025?
- Modify the input of these terms in `step.F90`.

ORCA025 has run without error for one month. SUCCESS.


### 10th July
... almost. ORCA025 crashed.

Reverting back to stretching calculations at every timestep, but using now temp and sal fields.

ORCA025 blows up with large velocities.


### 11th July
Have written in daily stretching calculations.
- Revision 16286

Model appears to break at first timestep...
ASCII output not working... Lets check idealised configuration.
	- Works for idealised... but `nmlnqg` is giving rubbish values. 

`IDEAL` appears to be running fine with daily stretching calculations ... Model blew up with large velocities at i, j, k = 214, 4, 25
- Computing at every timestep seems to keep the model stable

Try quarter daily stretching calculations in `u-cx856`.


### 19th July
Examining the ORCA025 output shows process discontinuities in vorticity and divergence gradients. 
- Compute these terms differently?
- Are the scale factors incorrect? Shouldn't gradients in `ji` use `e1v`? and vice versa.

In 2D Leith scheme only, have updated scale factors in vorticity, and added diagnostics. And in QG Leith scheme.


### 31st July
Failed to include name of subroutine at end of stretching calculation e.g. 
`END SUBROUTINE ldf_dyn_str`


### 1st August
Adding in print statements like,
`print *,  ‘dwzmagsq at’, ji, jj, jk, ‘is’, dwzmagsq(ji,jj,jk)`
which print to `testing.output`.

Print out only when `jk==22` to avoid large ascii output.

Messed up reading in of `zstlim` and the values do indeed match up with the print statements.


### 17th August
Possible cause for model blowing up may be due to use of `now` timestep variables. 
- Changing to use of before, means modifying `ldfdyn.f90` and `step.f90`. Doing this only in QG Leith routine to start.

New NEMO revision is 16324.


### 23rd August
ORCA025 crashed with before. But this seems the most consistent.

Adding in the stability criterion for QG Leith:
- Am = deltamin^2 / (8 deltaT).
- min horizontal resolution in orca025 is 7 km https://archimer.ifremer.fr/doc/2006/publication-3514.pdf
- Using normal timestep and this resolution gives a max viscous coefficient of 3000 m^2^/s, significantly lower than values being initially calculated.



### 19th September
We now begin to implement Leith as GM/Redi.

- In NEMO, GM coefficients are computed in `OCE/LDF/ldftra.f90`.
- There are two coefficient specified, one on each of u- and v- points.
- This can be achieved by averaging `ahmt` onto each point.

- Might need to re-order the calling of `ldfdyn.f90` in `step.f90`, so viscosity coefficient is called first then can be read into `ldftra.f90`.

To incorporate QG Leith, the easiest way may be to add IF conditions into `ldftra` to choose QG Leith when a namelist set TRUE in section `namtra_eiv`?


### 20th September
- When choosing Leith as GM, no need to compute the coefficient in routine `ldf_eiv` since it can be computed in `ldf_tra`. 
- The coefficients `aei` are assigned in `ldf_tra` anyway.

See commits on GitHub for further implementation details.

NEMO has compiled and IDEAL has run for 1 year with GM/Redi turned on.


### 28th September
Backtracking a little due to the QG Leith starry night sky.

We have attempted a few 'fixes' in `IDEAL`, but to no avail. Instead, we turn our attention to the numerical implementation of the stretching term, since this is where the starry night sky is appearing.

Instead of computing an x and y component of stretching, we compute one value. So following MITgcm implementation.


### 29th September
In an attempt to solve the starry night sky, we have implemented QG Leith in a similar fashion to MITgcm. A single stretching term computed on the T-point and then its gradient taken and averaged onto the vorticity gradient points, before the magnitudes are put onto T-point.

Preliminary output of IDEAL after 1 year from rest is not promising, where single grid point large stretching values are causing large gradients and grid box like viscosity values. 

We wonder if this is less to do with the implementation, perhaps more to do with spin up? Still the same.

Revert back to original `ldfdyn.F90`, try from spin up though.






