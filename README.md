----------------------------------------------------------------------
### Copyright (C) 2016 by Paolo Raiteri (Curtin University)

email: p.raiteri@curtin.edu.au                                        
                                                               
This program is free software; you can redistribute it and/or  modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.         
                                                               
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.                   
                                                               
For the full text of the GNU General Public License, write to: Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              
                                                               
The GNU GPL can also be found at http://www.gnu.org 
                                                               
No claim is made that this program is free from errors and no liability will be accepted for any loss or damage that may result. The user is responsible for checking the validity of their results. 

----------------------------------------------------------------------

# Reference

This program can be used to create disolcations pairs (Frank-Read dislocations) in molecular crystals as described in the paper:

Structure and dynamics of screw dislocations in even n-alkane crystals"

Authored by:
Isabel Olson,<sup>1</sup> Alexander G. Shtukenberg,<sup>1</sup> Gagik Hakobyan,<sup>1</sup> Andrew Rohl,<sup>2</sup> Paolo Raiteri,<sup>2</sup> Michael D. Ward,<sup>1</sup> Bart Kahr<sup>1,3</sup>

<sup>1</sup>Department of Chemistry and Molecular Design Institute, New York University, New York City, NY, 10003, USA.
<sup>2</sup>Curtin Institute for Computation and Department of Chemistry, Curtin University, P.O. Box U1987, Perth, Western Australia, 6845, Australia. 
<sup>3</sup>Graduate School of Advanced Science and Engineering (TWIns), Waseda University, Tokyo, Japan

published in XXX

----------------------------------------------------------------------

# Compilation

The code can be compiled by using the provided `Makefile` in the `src` directory.
Typical compilation options for *gfortran* and *ifort* have been provided, but you may have to edit the `Makefile` to suit your system.
If you have the *gfortran* compiler installed you can move into the `src` directory and simply type
```bash
make gf
```

to compile the code.
There is also a *debug* option that may help you trace back any issues with the code:
```bash
make gf-dbg
```

The executable `disloc.x` will be placed in the `bin` folder and after a succesful compilation all the object files will be moved to the `obj` folder, both located at the same level as `src`

An example of usage is provided in the `test` folder and it can run by typing
```bash
../bin/disloc.x dislocation.inp
```


----------------------------------------------------------------------

# Features


* The code requires an input file for the commands and a coordinate file in PDB format, see `test/dislocation.inp`
* The molecules are displaced along the burger vector using this function
$d=\pm \frac{1-\bigg(\displaystyle\frac{x}{\sigma}\bigg)^n}{1-\bigg(\displaystyle\frac{x}{\sigma}\bigg)^m}$
where $\sigma$, $n$ and $m$ are input parameters and $x$ is the distance between the centre of the disclocation (specified in the input) and the centre of mass of the molecule. Positive and negative displacements are applied to molecules sitting above and below the dislocation line, respectively.
* The handedness of the dislocation can be changed by inverting the direction of the dislocation line
