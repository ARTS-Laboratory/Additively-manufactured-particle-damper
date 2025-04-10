# simulation quick-start

## running YADE
Yade simulations are run via the terminal generally via:
> yadedaily [YADE code]

For these simulations specifically, use:
> yadedaily -j[# cores] [YADE code] [input file] [input file row]

the "-j[# cores] can be left out if no parallelization is wanted

If running with an indent, add "indent" to the end. Otherwise leave out of command.
> yadedaily -j[# cores] [YADE code] [input file] [input file row] indent

## setting input parameters
Input parameters are defined in a csv file. Currently utennessee.csv is being used, but any file name will work so long as the same headers are used.

## particle settlement (first simulation)
This generates the box and particles based on the input parameters, then lets the particles settle. 
> yadedaily -j[# cores] ball_gen_8_and_15.py [input file] [input file row] indent

Output files will be generated automatically over time to be used in the pocket agitation simulation.

## pocket agitation (second simulation)
This agitates the pocket using the same input parameters as the particle settlement simulation.
> yadedaily -j[# cores] yade_sdoftrial9.py [input file] [input file row] indent

Outpute files will be generated automatically over time.

## extra notes
A file called spec_vel.csv was used for research predating this repository. It contains potentially useful insights into the simulation's development.
A file called Run Line.txt exists and is used to store commands to be easily copy and pasted into the terminal.
The only files necessary prior to simulation are ball_gen_8_and_15.py, yade_sdoftrial9.py, common_functions.py, and an input csv.