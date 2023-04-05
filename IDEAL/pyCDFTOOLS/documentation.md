## Documentation of additional pyCDFTOOLS
This follows the work of J. Mak and CDFTOOLS/{}.f90

### cdftransport.py
### 28th March 2023
- This function computes the zonal transport of an arbitrary section.
- This has only been tested in a channel model setup.
- Want to add additional options to compute heat transports, with an option to calculate vertical heat transport.

### 29th March 2023
- Adding additional arguments into the function to enable the computation of heat transport, horizontally and vertically.
- Seems to be spitting out not so unreasonable numbers for vertical heat transport i.e. on the order of PW.
- Need to speed the loops up, perhaps using @jit method J.Mak used...

### 31st March 2023
- Making input arguments lists to enable multiple files and variables to be input. This comes from stackoverflow question. Must thank for responses.
- Renamed working copy of cdftransport that uses `**kwargs` `cdftransport.py.orig`.

Need to make updated `cdftransport.py` work.

### 5th April 2023
- Updated elif statement on line 80.
- Updated function description.

