## QG Leith Implementation
Should have started this a long time ago... 

Currently the routine is correctly called using the namelist parameters. The routine leans heavily on Scott Bachman's implementation in MITgcm.

See also the `documentation.md` file in `IDEAL` repo [here](https://github.com/thomaswilder/nemo-IDEAL) for further details on model setup to test the Leith schemes.

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
