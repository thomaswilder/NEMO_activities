# Objective
Documenting the implementation of the Leith schemes in ORCA025. We intend to run and modify the GOSI9p8.0 ORCA025 suite on Rose.

The idea behind this documentation is for me to document everything I do, from mistakes to useful procedures. Maybe this will be useful to someone else setting up a Rose suite on Monsoon2

### Useful resources
- https://code.metoffice.gov.uk/trac/moci/wiki/tips_CRgeneral
- Use `u-cn052(M)` on Monsoon2.
- Check [here](https://code.metoffice.gov.uk/trac/home/wiki/ProjectList) for svn url's.


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
