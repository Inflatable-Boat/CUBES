CUBES: Cubes Umbrella Biasing Entropic Simulation
=================================================

Usage of this repository
========================

In the main directory, compile (on Linux) with `gcc -Wall -O3 code/v21_sim.c -lm -o v21_sim.exe`. On Windows, `-lm` is not necessary.  
`-Wall` for all warnings, `-O3` for maximum optimization, `-lm` to link math library `math.h`.  
Note the program also needs the files `math_3d.h` and `mt19937.h` (also in this repository).

The output program `v21_sim.exe` will take 7 arguments:  
  1. (r or read / c or create)
  2. (readfile / Number of cubes per dimension)
  3. (output_folder / packing_fraction*)
  4. mc_steps, the number of Monte Carlo sweeps (= steps per particle)
  5. output_steps, every _output\_steps_ the program will make a snapshot of the system
  6. BetaP / NVT, the (dimensionless) pressure or do an NVT ensemble
  7. Phi, the slant angle

*= if there is overlap the system will expand (and thus reduce packing fraction).

The first three options are all left or all right, i.e. either 
  1. _create_ a system with _this number of cubes_, starting at _this packing fraction_, OR  
  2. _read_ a system from datafolder/<_readfile_> and save subsequent moves to datafolder/<_output_folder_>.

---
## Output

The program outputs some debug data at the start and then every _output\_steps_ it outputs the step number, the volume, step acceptance ratios (of move, rotation, and volume moves), and then the corresponding maximum stepsize of each kind of step (i.e. max distance per trial move, max rotation in radians per trial move, max change in volume per trial move).

It also outputs the density of the system every _output\_steps_ steps to e.g. `densities/v21_12pf0.50p01.0a1.27`.
(meaning 12 cubes per dimension, starting packing fraction of 0.50, pressure of 1.0, slant angle of 1.27 radians (~72 degrees))  
It then also outputs the radial distribution function (200 bins, from 0 to 5) of the system at this moment to a new line in `datafolder/v21_12pf0.50p01.0a1.27/g.txt`. Note: the first number is the number of particles counted.

Then it lastly also outputs snapshots of the system every _output\_steps_ steps into `datafolder/v21_12pf0.50p01.0a1.27/coords_stepXXXXXXX.poly`, where XXXXXXX is the step number.

---
## Examples

Example 1: `./v21_sim.exe c 4 0.80 10000 100 nvt 1.57079632679`  
**c**reates a system of **4** cubes per dimension,  
start at a packing fraction of **0.80**,  
runs for **10000** steps,  
make a snapshot every **100** steps,  
do an **NVT** simulation,  
set slant angle of cubes to **1.57079632679** radians (90 degrees).

Example 2: `./v21_sim.exe r datafolder/v21_04pf0.80p-1.0a1.57nvt/coords_step0010000.poly vtemp 10000 100 nvt 1.57079632679`  
**r**eads file **datafolder/v21_04pf0.80p-1.0a1.57nvt/coords_step0010000.poly**, runs for **10000** steps, saves a snapshot every **100** steps into datafolder/vtemp.

---
## Formatting of .poly
The format for the coordinate files is as follows.
For a selection of N slanted cubes of angle PHI and edge length L in a box described
by the matrix B (that is to say, the box is a parallelepiped spanned by the three
vectors that are the columns of the matrix B), their positions given by the vector R
and their orientation described by the rotation matrix M, the format is:

N  
0 0 0  
B_00 B_01 B_02  
B_10 B_11 B_12  
B_20 B_21 B_22  
R_0 R_1 R_2 L M_00 M_01 M_02 M_10 M_11 M_12 M_20 M_21 M_22 10 PHI COLOR  
 
Note that the rotation matrix, as its name implies, only describes a rotation of the 
particle from a predefined orientation. This predefined orientation is one where the
slanted cube's shape is given by:
v0 = ( L, 0, 0 )
v1 = ( 0, L, 0 )
v2 = ( L cos(PHI), 0, L sin(PHI) )

## Example:

The following would represent two slanted cubes, one placed in the origin, of edge length 1, unrotated from its default orientation, with a slant angle of 1.4 radians, and smaller cube placed at (1,1,1).

2  
0.000000 0.000000 0.000000  
2.371260 0.000000 0.000000  
0.000000 2.371260 0.000000  
1.025854 0.000000 2.137872  
0.0 0.0 0.0 1.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 10 1.4  
1.0 1.0 1.0 0.500000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 10 1.570796  

--- 
--- 
--- 
--- 
--- 
TODO
====

Heading
=======

## Sub-heading

Paragraphs are separated
by a blank line.

Two spaces at the end of a line  
produces a line break.

Text attributes _italic_, 
**bold**, `monospace`.

Horizontal rule:

---

Bullet list:

  * apples
  * oranges
  * pears

Numbered list:

  1. wash
  2. rinse
  3. repeat

A [link][example].

  [example]: http://example.com

![Image](clus2.png "icon")

> Markdown uses email-style > characters for blockquoting.