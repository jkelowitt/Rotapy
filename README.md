# Rotapy

![Logo](https://i.imgur.com/59TSCMn.png)

### Purpose

Generate Gaussian 09 input files for the rotamers of an input compound.

Distance to the axis of rotation remains constant throughout the rotation.

### Usage

This methyl group will be used as the example. Our goal is to make 6 rotomers with the Hydrogens shifted by 60 degrees
from each other.

![Image](https://i.imgur.com/e2AES96.png)

1a. If running the python directly, place the .log files you wish to 'rotamate' into the root directory of Rotapy.py,
and run Rotapy.py with python version >=3.6

1b. If running the executable, place Rotapy.exe in the same location as the .log files you wish to rotamate, and run
Rotapy.exe. When the terminal appears, you may need to press enter to get it started.

2. The program will ask for the directory of a .log Gaussian geometry optimization output file.


3. The program asks for an Anchor Atom. This is atom `1` as shown above.


4. The program asks for a Center atom. This is atom `2` as shown above.
    - The bond between the anchor and center atoms form the axis of rotation.
    - The center atom and the anchor atom must be different atoms.
    - The atoms attached to the center atom, `3, 4, 5`, will be rotated from their original position about the axis of
      rotation.


6. The program asks for an angle to scan through. To get rotomers with 60° increments, type `60`.
    - This will make rotomers representing 0°, 60°, 120°, 180°, 240°, and 300°.

    - Since 0° = 360°, the final rotamer is skipped.


7. The program asks if we want to add more rotations. Type `n`.
    - If we desired to generate rotamers of our rotamers, we would type `y` instead, repeating steps 2-5.
    - Two rotations of 60° results in `6 * 6 = 36` final rotamers.
    - Three rotations of 60°, plus one rotation of 10° results in `6 * 6 * 6 * 36 = 7,776` final rotamers.
    - All possible combinations of rotamers are made.


9. The program will ask if you want to save the rotamers to .com files, for processing, and will allow you to change the
   parameters for the Gaussian analysis. The editable parameters are:
    - title (Default: The name of the molecule as described below)
    - charge (Default: 0)
    - multiplicity (Default: 1)
    - job (Default: Opt Freq)
    - theory (Default: B3LPY)
    - basis_set (Default: 6-311G(2df,2p))
    - cores (Default: 8)
    - memory (Default: 20gb)
    - linda (Default: 1)


10. The program will ask if you want to save images of the rotamers to .png files.


11. The program gives one final warning, then performs the saving of .com and .png files to designated directories.

8. Once `n` is selected for step 6, the program performs, and prints out the name of each molecule in a specific format.
   - `[compound_name]_a[anchor atom number]c[center atom number]d[current value of angle rotation]`
   - For each rotation, another __a#-c#-##deg is appended to the rotamer name. For example:
      - `terpineol4_a7c17d60_a7c27d170_a17c18d180`
      - The compound Terpineol4 has:
         - The atoms on the 17th atom rotated by 60° about the axis formed between the 7th and 17th atom.
         - The atoms on the 27th atom rotated by 170° about the axis formed between the 7th and 27th atom.
         - The atoms on the 18th atom rotated by 180° about the axis formed between the 17th and 18th atom.

## Work remaining

- Check that no bonds are formed during the rotation, and skip over the rotamer where they form.
- A GUI. ["artist's" representation](https://puu.sh/HJ8K5/73b7ca6259.jpg)
- Allow for individual rotations, rather than a scanned rotation.

### References used in code

- https://github.com/tmpchem/computational_chemistry
- https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
