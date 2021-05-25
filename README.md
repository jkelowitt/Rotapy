# Rotapy

### Purpose

Generate Gaussian 09 input files for the rotamers of an input compound.

Distance to the axis of rotation remains constant throughout the rotation.

### Usage

This methyl group will be used as the example. Our goal is to make 6 rotomers with the Hydrogens shifted by 60 degrees
from each other.

![Image](https://i.imgur.com/e2AES96.png)

1. Run `main.py` using Python 3.6 or higher.

2. The program will ask for the directory of a .log Gaussian geometry optimization output file.

3. The program asks for an Anchor Atom. This is atom `1` as shown above.


4. The program asks for a Center atom. This is atom `2` as shown above.
    - The bond between the anchor and center atoms form the axis of rotation.
    - The center atom and the anchor atom must be different atoms.


5. The program asks for the atoms to be rotated. These atoms are the Hydrogens. They are input as `3,4,5`, with no
   spaces between the numbers.
    - The atoms to be rotated don't necessarily need to be bonded to the center or anchor atoms.


6. The program asks for an angle to scan through. To get rotomers with 60° increments, type `60`.
    - This will make rotomers representing 0°, 60°, 120°, 180°, 240°, and 300°.

    - Since 0° = 360°, the final rotamer is skipped.


7. The program asks if we want to add more rotations. Type `n`.
    - If we desired to generate rotamers of our rotamers, we would type `y` instead, repeating steps 2-5.
    - Two rotations of 60° results in `6 * 6 = 36` final rotamers.
    - Three rotations of 60°, plus one rotation of 10° results in `6 * 6 * 6 * 36 = 7,776` final rotamers.
    - All possible combinations of rotamers are made.


8. Once `n` is selected for step 6, the program performs, and prints out the name of each molecule in a specific format.
    - `[compound_name]__a[anchor atom number]-c[center atom number]-[current value of angle rotation]deg`
    - For each rotation, another __a#-c#-##deg is appended to the rotamer name. For example:
        - `terpineol4__a7-c17-60deg__a7-c27-170deg__a17-c18-180deg`
        - The compound Terpineol4 has:
            - The atoms on the 17th atom rotated by 60° about the axis formed between the 7th and 17th atom.
            - The atoms on the 27th atom rotated by 170° about the axis formed between the 7th and 27th atom.
            - The atoms on the 18th atom rotated by 180° about the axis formed between the 17th and 18th atom.


9. The program will ask if you want to save the rotamers to .com files, for processing, and will allow you to change the
   values for the analysis.


10. The program will ask if you want to save images of the rotamers to .png files.


11. The program gives one final warning, then performs the saving of .com and .png files to designated directories.

## Work remaining

- Recursively search for atoms attached to the center atom, so the user doesn't have to individually add them.
- Check that no bonds are formed during the rotation, and skip over the rotamer where they form.
- A GUI. "[artist's representation](https://puu.sh/HJ8K5/73b7ca6259.jpg)"
- Allow for individual rotations, rather than a scanned rotation.

### References used in code

- https://github.com/tmpchem/computational_chemistry
- https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
