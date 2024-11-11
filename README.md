**<h1>Simulated Annealing Engine for Slicing Floorplan**

**TEAM: IDEAL IDLE</h1>**


A Python-based tool called Simulated Annealing Engine for Slicing Floorplan
is intended to utilize a SA algorithm for providing a sub-optimal positioning
of any given blocks (HARD & SOFT MACROS) for any given area of floorplan.
It performs the following functions

1. a) **Extraction:** of block values from the .blocks (HARD & SOFT MACROS) GSRC benchmark file.

   b) **Computation:** Area of Floorplan and the near-optimal coordinates of the blocks.

	
**<h3>DEPENDENCIES:</h3>**

-   Python 3.7 or higher

-   Libraries:

    -   '**argparse**' for parsing command-line arguments.

    -   '**math**' for computing math operations.
	
	-   '**random**' for generating random numbers.
	
	-   '**functools**' for caching information.

**<h3>USAGE:</h3>**

- Bash Terminal Commands Supported for a given mode:

**(1) Instructions to run floorplanner.py**

**To compute floorplan from a .blocks file:**

>> python3.7 floorplanner.py -input HARD/n10.blocks -output n10.out (path to .blocks file)

[Change source and output name as per (.blocks) file requirement]

Function: Provides a .out text file with processed coordinates of blocks with final and unused
computed floorplanning area fetched from the SOFT/HARD .blocks File.

**<h3>EXAMPLE(S):</h3>**

**(a)To analyze a benchmark file named n10.blocks and compute floorplan details:**

>> python3.7 floorplanner.py -input HARD/n10.blocks -output n10.out

**<h3>OUTPUTS</h3>**

The generated output files after running any of the above Bash Commands are stored under the folder '**outputs**' in the parent directory.
