# Objective
Documenting the implementation of the Leith schemes in ORCA025. We intend to run and modify the GOSI9p8.0 ORCA025 suite on Rose.

The idea behind this documentation is for me to document everything I do, from mistakes to useful procedures. Maybe this will be useful to someone else setting up a Rose suite on Monsoon2

### Useful resources
- https://code.metoffice.gov.uk/trac/moci/wiki/tips_CRgeneral
- Use `u-cn052(M)` on Monsoon2.
- Check [here](https://code.metoffice.gov.uk/trac/home/wiki/ProjectList) for svn url's.
- [Getting started](https://code.metoffice.gov.uk/trac/home/wiki/FAQ).


### 21st April
- Checked out a local copy of `u-cn052`, which is called `u-cr756` in `rosie go`.
- Trying to run the suite using `rose suite-run --new`, but not working?
	- Some tasks in cycl GUI stuck on `submit-retrying`.
	- Looking at their log files e.g. for `fcm_make2_ocean` in `job-activity.log` shows `(xcslc0) 2023-04-21T13:52:09Z [STDERR] qsub: error: [PBSInvalidProject] 'jmmp' is not valid for collaboration trustzone on XCS`.
- Under `suite conf -> UM hosts -> PROJECT_USER` change `jmmp` to `ukesm`.
	- rose suite seems to be running now.

What svn url do i create a branch from to work on `u-cn052`? Is it here https://code.metoffice.gov.uk/svn/nemo/  ?

Errors in postprocessing of `u-cr756`.
- Seems to still be trying to use `jmmp` project when executing moose command.
- Shutting down suite for the day.


### 24th April
- Running u-cr756 again to reproduce error above.

Trying to locate GOSI9p8.0:
	- Located [here](https://code.metoffice.gov.uk/trac/nemo/log/NEMO/branches/UKMO/NEMO_4.0.4_GOSI9_package?rev=16215) on trac.
	- Located [here](https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/UKMO/NEMO_4.0.4_GOSI9_package/) on svn.
cr756o_1m_19760201_19760228_grid_T_197602-197602.nc
Going to play around checking out and making changes with the test repository provided [here](https://code.metoffice.gov.uk/trac/test).
	- Test ticket number 32.
	- `Committed revision 232.
[info] Created: https://code.metoffice.gov.uk/svn/test/test/branches/dev/thomaswilder/r231_wildert`

- Can find the source code location using `rosie go` and going to `fcm_make_ocean -> env -> NEMO and SI3 sources`.
- Need to make my own branch of the GO8 package.

Created a ticket on https://code.metoffice.gov.uk/trac/GO/ entitled 'Implementing Quasi-Geostrophic Leith Viscosity'. 
	- Ticket number 652.
	- Location: https://code.metoffice.gov.uk/trac/gmed/ticket/652#ticket

Checked out a branch of GO8 located here https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/dev/thomaswilder/r15557_NEMO_4.0.4_GO8_package_qgleith

To do:
1) Incorporate QG Leith code change into working copy of branch.
2) Commit changes in working copy to metoffice branch.
3) Test QG Leith.

First, why is postproc task not completing.
- In rosie go, `postproc -> Post Processing - common settings -> Moose Archiving`, change `moo-project` to `project-ukesm`.
- Try running for 2 months.
- Changed `base_component` to `5 days` under `postproc -> NEMO -> Diagnostics -> Meaning`.

### 25th April
Saved `moo-project` to `project-ukesm` in `postproc`.

Realised I created a nemo branch from the trunk. Deleted branch, and creating a new one using the `--branch-of-branch` option in `fcm create-branch`.

Created again using the NEMO trunk!!! Why? Maybe its a sub tree so cannot create a branch from that?
Tried also creating a branch from trunk at an earlier revision.
- Queried this with Till by email.

ORCA025 ran successfully for one month and wrote data to file.
	- Data restart and grid files located in `cycl_run/u-cr756/share/data/History_Data/NEMOhist`
	- Checking if data written to mass with `moo ls moose:/crum/u-cr756`. Files are `ida.file`, `oda.file`, and `onm.nc.file`...?
	- `onm.nc.file` contains all nemo grid files.

Installed moose on JASMIN.
	- Log in like `ssh -X twilder@mass-cli.jasmin.ac.uk`
	- Can also view files of `u-cr756` from mass VM on JASMIN in the same way as is done on Monsoon2.

How do I add diagnostics to the ORCA025 suite?
	- Can fine the `file_def_nemo-oce.xml` in `/roses/u-cr756/app/nemo/file/` directory.	

### 27th April
Can we restart the rose suite by changing the pickup/start date?
	- Lets also turn off compiling drivers and ocean.
Seems to be running.


### 28th April
Picked up from month 1 and ran for a further month successfully.

Try and modify the diagnostics.
	- Potentially do this by selecting custom diagnostics and including the `file_def_nemo-oce.xml` in `/file` directory of rose suite.
	- In `file_def_nemo-oce.xml`, why is the first `file id="file8"`?
	- So, `file_def` files are in `roses` directory, but `field_def` is not, so how do you add diagnostics to the rose suite? Add your own file and point the suite to that file? Put in same `file` directory as the `file_def`.

To Do:
1) Modify `file_def_nemo-{oce,ice}.xml` files.
2) Will the rose suite run?


### 2nd May
To Do:
1) See email about duplex data in mass.


### 18th May
To Do:
1. Checkout NEMO branch,
2. Add modified `ldfdyn.F90` and `step.F90` into NEMO source code,
3. Update `field_def_nemo-oce.xml` and `file_def`,

Found a possible workaround for creating a nemo branch [here](https://code.metoffice.gov.uk/trac/GO/wiki/GODocumentation/GO8.0/NemoUpgradeInstructions#no1)

Possibly successful:
- Made a copy of the branch `NEMO_4.0.4_GO8_package` into `dev/thomaswilder/...`. This is found [here](https://code.metoffice.gov.uk/trac/nemo/browser/NEMO/branches/dev/thomaswilder/NEMO_4.0.4_GO8_package_QGLeith?rev=16231).
- Need to check if the Rose suite will run with a link to the dev branch first before making any changes.
- `fcm bc` doesn't do what I want it to do.
- Have updated CMS post with solution.


### 24th May
Working on:
	- Modifying source code,
	- Modifying field and file def files.

The `file_def` and `field_def` files are extracted from the `METO_GO` cfgs on `code.metoffice` `GO8` branch.
	- Add own custom files to the `ocean_ice` folder in the roses suite.
	- Try putting BASIC diagnostic files in `app/file` and seeing if this works. Changing the file numbers. 
	- What does `sync_freq` mean in the xml fles?

Model has run and data has outputted, and looks good. Have turned off postproc and housekeeping for now.
	- Run postproc and housekeeping.
	- Added `eken` diagnostic to check if this works. Picking up from restart file and running only model for one month. This works.

To Do:
1. Modify `ldfdyn` for mixed layer depth, use the variable `nmln` from original source. Test this first in the idealised config.
2. Add in modified src files and commit these changes to the branch on MO,
3. Modify field def file to include QG diagnostics.


### 1st June
From ESMGA, consider adding in ice cavities to the GOSI9 config, and diagnose the impact of eddy parameterisations on them.


### 5th June
- Updated [ticket](https://code.metoffice.gov.uk/trac/gmed/ticket/652) with branch information.
- Adding in modified src files to the branch on MO ... done. Revision now updated to 16236. Need to account for this in rose suite.
- Lets see if QG Leith is an option in the rose suite ...? No. Try adding the options manually to `app/nemo/rose-app.conf`. Error in rose suite saying `value 34 not in allowed values`. So how how is it added to the choices... reference? But it might just work. Lets try it...
- Adding `ahmt_3d` diagnostic to `file_def_nemo-oce.xml`. Uncertain where `field_def` goes, maybe in same directory as `field_def`?

Does it compile? YES.

Now run ocean and postproc components. Fail on `ocean_ice`,
- `> Error [CObjectFactory::GetObject(const StdString & id)] : In file '/var/spool/jtmp/6413365.xcs00.XplTEw/xios/src/object_factory_impl.hpp', line 78 -> [ id = rro2, U = field ] object was not found.`

Lets try 2D Leith. Looks to be running fine... What do the output fields look like? Check JASMIN. Looks great.

Error in QG Leith could be due to diagnostics, `rro2` is the first `iom_put...`. Need to modify `field_def` accordingly.
- Uploading modified `field_def...` from `NEMO_activities/ORCA025` to `ocean_ice` app. 

MOOSE data:
- `ida.file` are restart files.
- `onm.nc.file` are ocean diagnostic files.

Testing out QG Leith with `field_def...` added to `ocean_ice` app ... This did not work.
- Maybe change source for file in rose suite. This works, but for some reason the `field_def` file from `cfgs/SHARED` on MO repo is not the same one extracted in `cylc-run`.
- Copied one from cylc directory, modified it, and copied in to `ocean_ice` app.
- Seems to run for a while but eventually fail with error,
```
Application 196964917 is crashing. ATP analysis proceeding...
Rank 320 [Mon Jun  5 15:40:06 2023] [c7-2c0s8n3] application called MPI_Abort(MPI_COMM_WORLD, 60) - process 320

ATP Stack walkback for Rank 320 starting:
  _start@start.S:113
  __libc_start_main@libc-start.c:242
  main@nemo.f90:18
  main@nemo.f90:18
  nemo_gcm$nemogcm_@nemogcm.F90:167
  stp$step_@step.F90:261
  stp_ctl$stpctl_@stpctl.F90:45
  ctl_stop$lib_mpp_@lib_mpp.F90:1722
  mppstop$lib_mpp_@lib_mpp.F90:1299
  mpi_abort__@0x56024bc
  MPI_Abort@0x55f97ed
  MPID_Abort@0x5635661
  abort@abort.c:92
  raise@pt-raise.c:42
ATP Stack walkback for Rank 320 done
Process died with signal 6: 'Aborted'
Forcing core dumps of ranks 320, 108, 144, 185, 115, 156, 180, 157, 184, 290, 301, 0, 344, 348, 346, 349, 216, 252, 333, 222
View application merged backtrace tree with: stat-view atpMergedBT.dot
You may need to: module load stat

_pmiu_daemon(SIGCHLD): [NID 06768] [c7-2c0s12n0] [Mon Jun  5 15:41:51 2023] PE RANK 342 exit signal Killed
[NID 06768] 2023-06-05 15:41:51 Apid 196964917: initiated application termination
[FAIL] run_model # return-code=137
2023-06-05T15:42:12Z CRITICAL - failed/EXIT
```

