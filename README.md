# Reference for triple.c
As used in https://doi.org/10.1073/pnas.2410109121:
Ding, E. A., Yokokura, T. J., Wang, R., Kumar, S. "Dissecting neurofilament tail sequence-phosphorylation-structure relationships with multicomponent reconstituted protein brushes." Proc. Natl. Acad. Sci. 121. 2024.

Code written by: Yokokura, T. J., Duan, C., Wang, R. UC Berkeley, Department of Chemical and Biomolecular Engineering. Primary contact: ruiwang325@berkeley.edu

## Files
1. `triple.c`: Source
2. `graft1.txt`: General input
3. `ac.txt`: Protein input

## Compilation and run
```
gcc triple.c -lm
./a.out
```

## Dependencies
1. Standard libraries (`stdio.h, stdlib.h, math.h, string.h`)
2. All files as listed above with `ac.txt` able to be edited as referenced in Line 2 of `graft1.txt`.
3. If initializing from previous output, `Win` as referenced in Line 12 of `graft1.txt`.

## Parameters
*Note: Descriptions for each line can be found in the comments of triple.c. The following is an abbreviated guide to run the code and does not offer solutions for debugging or suggestions for convergence, which both depend on the system at hand.* 

`graft1.txt` (assuming ternary mixture is used. If not, take out Line 5 for binary and Lines 4 & 5 for single-component): 
* Line 1: Initialization options
    * `-1`: Read fields from `Win` (Line 12)
    * `0`: Start from homogeneous polymer density of 1st value of Line 19
* Line 2: `filepath` for protein blocks
* Lines 3-5: `Grafting density` of each protein (up to 3)
* Line 6: box `size` [nm]
* Line 7: Discretization `number` along contour of Gaussian chain
* Line 9: Bulk salt `concentrations` [M]
* Line 15: Mixing parameters for solution convergence
* Line 16: Anderson mixing parameters
* Line 16: Convergence thresholds
* Line 20: Kuhn segment `[length` and `volume]` for each protein

`ac.txt` (filename determined by Line 2 of `graft1.txt`): 
* Line 1: proteins considered (`LMH` for ternary / `LM` for binary / `L` for single-component)
* Line 2: `Number` of blocks within protein
* Line 3-X: `[Start, end]` residue number for a given block, `charge density`, Flory--Huggins `Chi` parameter
* Repeat up to two more times for binary and ternary brushes