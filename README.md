# Rotapy

![Logo](https://i.imgur.com/59TSCMn.png)

### Purpose

Generate Gaussian 09 input files for the rotamers of an input compound.

Distance to the axis of rotation remains constant throughout the rotation.

### Usage

To walk through how to use Rotapy, we will walk through the rotation of the alcohol group of ethanol.

Start by opening Rotapy:

![Empty](https://puu.sh/HWOwv/23f6d24335.png)

In this area here, you can select a .com, .log, or .xyz file to analyze. This can be done by clicking Browse and finding
the file, or by copying and pasting the absolute file path into the input box:

![Selecting Import Molecule Browse](https://puu.sh/HWOQT/f9f0e7468a.png)

![With Import Molecule Text](https://puu.sh/HWOTI/5ef376a1d7.png)

With the file selected, click on the "Show Molecule" button to view the molecule. This view is 3d projected into 2d, and
doesn't have proper perspective. It is possible to rotate and view the molecule from different angles however.

![Figure View](https://puu.sh/HWOUS/e7e9847d98.png)

We will need to make note of two atoms. The anchor atom, assuming we are rotating the alcohol group, is atom 1. The
center atom, the atom which has all the substituents we want to rotate, is atom 5.

We can go back to the main rotapy window, and type in the number for the anchor atom in the anchor input, and the center
atom into the center input.

![Center and Anchor Added](https://puu.sh/HWOWy/0643e4fbbd.png)

If we wanted to make rotamers of the 15 degree rotations of the alcohol group, we type 15 into the Angle input.

![Angle Added](https://puu.sh/HWOXS/1e7d30d8e7.png)

Once all three numbers are selected, clicking the "Add" button will add the rotation to the queue.

![Add Selected](https://puu.sh/HWP0i/4d67d38399.png)
![Added Alcohol Rotation](https://puu.sh/HWP1p/de55d56ab1.png)

Notice that the Total Rotamers has increased to 24. This means that in the end, there will be 24 rotamers made. If we
add a rotation to the queue representing the methyl group being flipped 180 degrees, we can see the number of rotamers
double.

![Added Methyl Rotation](https://puu.sh/HWP8F/28593e9183.png)

If we wanted to remove a rotation from the queue, we can click on the rotation in the queue, then click the "Remove"
button.

![Selected Methyl Rotation](https://puu.sh/HWP8g/5c6dfa2e87.png)
![Removed Methyl Rotation](https://puu.sh/HWP3u/ea29d5f1b2.png)

Now that we have all the rotations that we want added, we can prepare for getting output. If we want to have .com files
as an output, a directory must be selected. If a directory is not in the input box, the .com files will not be
generated. The same goes for the image files. Images of the rotamers will not be generated if the input box is empty. A
path to the directory must be provided, either by pasting the absolute path to the directory, or by browsing to the
directory.

![Added Output Directories](https://puu.sh/HWPdN/76c8157022.png)

Since both the Com and Img outputs both have a directory, .com and .png files will be generated for the rotamers.

Now that we have the output locations, we can click the "Change Output Settings" button for more .com file output
options.

![Added Output Directories](https://puu.sh/HWPfg/3367a83cb4.png)

From here you can change the charge, multiplicity, job, level of theory, basis set, the number of cores, the amount of
memory, and the linda cores used. There is no validation for these inputs, so you may go through the process of rotating
a ton of files, and end up with invalid option errors in Gaussian. If this happens,
use [this program](https://github.com/jkelowitt/ConvertToCom) to generate new com files with different settings, without
having to perform the rotation calculations again.

There is also a setting called "Sequentially Name Files". If this setting is used, the files will be labelled, file_1,
file_2, etc. If this setting is not checked off, the explicit rotations will be added to the file name instead. In the
example of the methyl and alcohol group of ethanol being rotated, this would look like: ethanol_1a0c180d,
ethanol_1a0c180d_1a5c30d, etc. 1a0c180d means, anchor atom 1, center atom 0, rotated 180 degrees. If there are multiple
rotations in the same file name, they are listed. 1a0c180d_1a5c30d means that we rotate the methyl group 180 degrees,
then the alcohol group 30 degrees. This results in very long file names, which is why this is not the default setting.

The "Reset to Default" button will reset the values of the inputs to the settings shown above.

Clicking the "Save" button will close the window, but will save the options you have chosen. The saved options are lost
when Rotapy is closed, or the "Reset to Default" button is clicked.

With the options selected, click the "Perform Calculations" button to execute the rotation queue.

![Added Output Directories](https://puu.sh/HWPsS/bc48bc85fc.png)

The progress bar will show what step Rotapy is working on, and how much progress is left. Once it has finished, Rotapy
will pop up a message saying it is done.

![Added Output Directories](https://puu.sh/HWPui/51532aa479.png)

If given a directory, the com and image files will be contained within, named with the selected naming scheme, and
formatted with the selected options.

## Work remaining for me

- Allow for individual rotations, rather than a scanned rotation.

## References used in code

- https://github.com/tmpchem/computational_chemistry
- https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
