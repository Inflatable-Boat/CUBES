TODO
===

The format for the coordinate files that this opengl code visualizes is as follows.
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

For example (this is test.poly):

1
0.000000 0.000000 0.000000
2.371260 0.000000 0.000000
0.000000 2.371260 0.000000
1.025854 0.000000 2.137872
0.0 0.0 0.0 1.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 10 1.4 


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

![Image](Image_icon.png "icon")

> Markdown uses email-style > characters for blockquoting.