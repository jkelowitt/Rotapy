"""
@Author: Jackson K Elowitt
@Start Date: May 14, 2021
@Contact: jkelowitt@protonmail.com
@Site: github.com/jkelowitt/Rotapy

The end goal of this script is to be able to take in a Gaussian09 .log file,
and allow the user to rotate specific elements of the contained molecule, and
generate Gaussian09 input files for optimization on a supercomputer. Ideally,
multiple rotations may be performed simultaneously, such as in a nested for loop.

Main Changes remaining:
    - Theres a lot of faffing about switching between atom numbers and the atom itself.
        There's got to be a better way to handle this.
        Try looking up mathematica graph based code to find a better way to code Molecule
    - Check that no bonds break or add during the rotation.

"""

import multiprocessing
import sys
from functools import partial
from glob import glob
from multiprocessing import Pool

from numpy import arange
from tqdm import tqdm

from classes import Molecule, Atom
from functions import (show_structure,
                       save_structure,
                       bonded_atom_search,
                       center_on_atom,
                       rotate_point_around_vector,
                       verified_input)
from parsing import (make_choice_list,
                     yes_no,
                     parse_opt_geom_from_log,
                     make_choice_dict,
                     write_job_to_com)

# Prevents errors with the executable
multiprocessing.freeze_support()


def main():
    """
    Main program function

    Asks the user for the file they want to make rotamers out of.
    Asks the user for the anchor atom
    Asks the user for the center atom
    Asks the user for the rotation angle step size
    Asks the user if they'd like to add more rotations
    Asks the user if they want to save generated .com files
    Asks the user if they want to save generated .png files
    Performs the rotations, and saves any specified files.

    """

    # Contains info on how to perform all desired rotations
    # Sublists have the format:
    # [center_atom, ancr_atom, [rotatees], angle increments]
    rotation_queue = []
    rotation_count = 1

    # Get all the log files in the current directory and all subdirectories.
    files = glob("*.log")

    # Check that there are log files to be found.
    if not files:
        input("No log files found in the current directory or lower.")
        sys.exit()

    choice = make_choice_list(files)
    name_xyz = parse_opt_geom_from_log(choice)
    atoms = [Atom(a[0], (a[1], a[2], a[3])) for a in name_xyz]
    molecule_name = input("What is the name of this compound: ")

    base_compound = Molecule(molecule_name, atoms)

    print("You're going to be asked for an anchor and a center atom.")
    print("For the example of rotating an alcohol:")
    print("    C ----- O ----- H    ")
    print("    ^       ^       ^    ")
    print("    Anchor  Center  Will be rotated\n")

    while True:
        # Ask before showing plot
        if yes_no("View numbered structure: "):
            show_structure(base_compound)

        ancr_atom_num = verified_input("Which atom is the anchor atom (ex. 7): ", int)
        center_atom_num = verified_input("Which atom is the center atom (ex. 27): ", int)

        base_compound = center_on_atom(base_compound, center_atom_num)
        ancr_atom = base_compound.atoms[ancr_atom_num]
        center_atom = base_compound.atoms[center_atom_num]

        # All things attached to the center atom will be rotated.
        # The first atom is the center atom. Remove it to prevent numbering confusion.
        rotate_atoms = bonded_atom_search(base_compound, start=center_atom, wall=[ancr_atom])[1:]

        # Remove duplicates (and order btw)
        rotate_atoms = list(set(rotate_atoms))

        # Convert from Atom to number again.
        rotate_nums = [base_compound.atoms.index(atom) for atom in rotate_atoms]

        # Get rotation amount
        print(
            f"\nThe {len(rotate_nums)} atom(s) attached to the center atom will be scanned through their rotation."
        )
        angle = verified_input("What step size should the scan perform, in degrees (ex. 45deg -> 8 rotamers): ", float)
        degrees = arange(angle, 360, angle)

        rotation_queue.append(
            {
                "center": center_atom_num,
                "ancr": ancr_atom_num,
                "rotatees": rotate_nums,
                "angles": degrees,
            }
        )

        # Increase total rotation count
        rotation_count *= len(degrees) + 1

        print(
            f"\nThe rotation of {len(rotate_nums)} atoms,\n"
            f"\tabout the {center_atom_num},{ancr_atom_num} atom bond axis,\n"
            f"\twith a total of {len(degrees) + 1}, {angle} degree rotations,\n"
            f"\thas been added to the queue of rotations.\n"
            f"\tThere are now {rotation_count} rotations total.\n"
        )

        if not yes_no("Add more rotations"):
            break

    rotamers = [base_compound]

    # Options on what and where to save the rotamers.
    # Default save location
    com_output = ""
    image_output = ""

    # Get saving preferences
    save_com_files = yes_no("\nSave the rotations to .com files")
    save_images = yes_no("Save images of the the rotations to .png files")

    # Get saving locations
    if save_com_files:
        com_output = input("\nWhat would you like to name the output directory for the com files: ")

        settings = {
            "charge": "0",
            "mul": "1",
            "job": "Opt Freq",
            "theory": "B3LPY",
            "basis set": "6-311G(2df,2p)",
            "cores": "8",
            "mem": "20gb",
            "linda": "1",
        }

        # Display default settings
        print("\nDefault Settings:")
        for item in settings:
            print(f"\t{item} = {settings[item]}")

        non_default = yes_no("\nUse the default settings")

        # Change default options if desired
        if not non_default:
            settings["Exit"] = ""
            while True:
                option = make_choice_dict(settings, prompt="\nWhich would you like to change (ex. job): ")

                if option in ("exit", "Exit"):
                    break

                value = input("\nWhat would you like to change it to (ex. opt): ")
                settings[option] = value
                print(f"Changed {option} to {value}\n")

    if save_images:
        image_output = input("\nWhat would you like to name the output directory for the image files: ")

    # Pause for preparation and alert the user to the download file count.
    input(f"\nPress enter to perform the {rotation_count} rotations and save the files.")

    # Perform rotation calculations

    erroneous_rotations = 0

    # With a progress bar
    with tqdm(total=rotation_count, desc="Performing rotation calculations", dynamic_ncols=True) as pbar:

        # Step through the rotation queue
        for rotation in rotation_queue:

            # Counter holds the number of molecules to have rotamers made from.
            # This is done to prevent performing rotations twice on the same rotamer.
            counter = len(rotamers)

            # Step through the molecule list
            for count in range(counter):

                # Step through the angles to be performed
                for turn in rotation["angles"]:
                    # Rotated atoms
                    rotated = []

                    # The molecule being rotated
                    originator = rotamers[count]

                    # Center the molecule on the center atom
                    new_rotamer = center_on_atom(originator, rotation["center"])

                    # Get the ancr atom to represent the axis of rotation
                    ancr_atom = new_rotamer.atoms[rotation["ancr"]]

                    # For every atom we want to rotate around this particular axis
                    for atom_num in rotation["rotatees"]:
                        # Get the atom to be rotated
                        atom = new_rotamer.atoms[atom_num]

                        # Find it's new position
                        new_pos = rotate_point_around_vector(atom.pos, ancr_atom.pos, turn)

                        # Add the new atom to the finished atoms list
                        rotated.append(Atom(atom.name, new_pos))

                    # Replace the un-rotated atoms with the rotated atoms
                    for old, new in zip(rotation["rotatees"], rotated):
                        new_rotamer.replace_atom(old, new)

                    new_rotamer.name += f"__a{rotation['ancr']}-c{rotation['center']}-{turn}deg"

                    # Due to how the rotation is performed, bonds will never break, only form.
                    # Bonds being added indicates that two atoms have gotten too close to each other.
                    # The compound is not a rotamer if new bonds are formed, so exclude these cases.

                    # Makes a list of the bond counts for each atom
                    originator_bond_count = [len(originator.bonds[b]) for b in originator.bonds]
                    new_bond_count = [len(new_rotamer.bonds[b]) for b in new_rotamer.bonds]

                    originator_bond_count.sort()
                    new_bond_count.sort()

                    # If the bond count lists are not exactly equal, a collision occurred
                    # Only a very complex scenario would really defeat this detection method.
                    # A complex problem == complex solution, thus, I procrastinate.
                    if new_bond_count == originator_bond_count:
                        rotamers.append(new_rotamer)
                    else:
                        bond_diff = sum(new_bond_count) - sum(originator_bond_count)
                        new_rotamer.name += f"_{bond_diff}ERR"
                        rotamers.append(new_rotamer)

                        erroneous_rotations += 1

                    pbar.update(1)

    # Alert the user to the presence of bad rotamers
    if erroneous_rotations:
        print(f"During the rotation calculations, {erroneous_rotations} rotamers formed or broke bonds. "
              f"They will be labelled with 'ERR'.")

    # Perform file saving
    if save_com_files:
        for molecule in tqdm(rotamers, desc="Saving com files", dynamic_ncols=True):
            write_job_to_com(molecule.atoms, title=molecule.name, output=com_output)

    if save_images:
        # args and kwargs can't be passed into the imap function, so we make a partial
        # function with the kwargs passed in here.
        saving = partial(save_structure, output=image_output)

        # Using pools results in ~4x speed performance boost when saving images.
        # 360 rotations in 00:02:45 -> 00:00:42
        with Pool() as pool:
            # Using imap, despite being slower than pool.map, so that we can have a progress bar.
            list(pool.imap(saving, tqdm(rotamers, desc="Saving images", dynamic_ncols=True)))


if __name__ == "__main__":
    main()
    input("Calculating and saving complete. Press enter to close. ")
