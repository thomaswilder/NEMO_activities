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
	- Potentially do this by selecting custom diagnostics and including the `file_def_nemo-oce.xml` in `/file' directory of rose suite.
	- In `file_def_nemo-oce.xml`, why is the first `file id="file8"`?
	- So, 	file_def` files are in `roses` directory, but `field_def` is not, so how do you add diagnostics to the rose suite? Add your own file and point the suite to that file? Put in same `file` directory as the `file_def`.

To Do:
1) Modify `file_def_nemo-{oce,ice}.xml` files.
2) Will the rose suite run?

