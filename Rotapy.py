"""
@Author: Jackson K Elowitt
@Start Date: May 14, 2021

The end goal of this script is to be able to take in a Gaussian09 .log file,
and allow the user to rotate specific elements of the contained molecule, and
generate Gaussian09 input files for optimization on a supercomputer. Ideally,
multiple rotations may be performed simultaneously, such as in a nested for loop.

Main Changes remaining:
    - Theres a lot of faffing about switching between atom numbers and the atom itself.
        There's got to be a better way to handle this.
    - Check that no bonds break or add during the rotation.

"""

import multiprocessing
import sys
from dataclasses import dataclass, field
from functools import partial
from glob import glob
from multiprocessing import Pool

import matplotlib
import matplotlib.pyplot as plt
import pyquaternion as pq
from numpy import arange, sqrt
from tqdm import tqdm

from data_dicts import cov_rads, bond_threshold
from parsing import (make_output_folder,
                     make_choice_list,
                     yes_no,
                     parse_opt_geom_from_log,
                     make_choice_dict,
                     write_job_to_com)

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

    # Ask before shwoing plot
    if yes_no("View numbered structure:"):
        plot_structure(base_compound)

    while True:
        ancr_atom_num = int(input("Which atom is the anchor atom (ex. 7): "))
        center_atom_num = int(input("Which atom is the center atom (ex. 27): "))

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
        angle = float(
            input(f"What step size should the scan perform, in degrees (ex. 45deg -> 8 rotamers): ")
        )
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
                option = make_choice_dict(
                    settings, prompt="\nWhich would you like to change (ex. job): "
                )

                if option in ("exit", "Exit"):
                    break

                value = input("\nWhat would you like to change it to (ex. opt): ")
                settings[option] = value
                print(f"Changed {option} to {value}\n")

    if save_images:
        image_output = input(
            "\nWhat would you like to name the output directory for the image files: "
        )

    # Pause for preparation and alert the user to the download file count.
    input(
        f"\nPress enter to perform the {rotation_count} rotations and save the files."
    )

    # Perform rotation calculations
    # TODO take this out and make it threaded. >50/sec is the goal
    with tqdm(total=rotation_count, desc="Performing rotation calculations") as pbar:
        for rotation in rotation_queue:  # Step through the rotation queue
            counter = len(rotamers)  # Counter holds the number of molecules to have rotamers made from.
            for count in range(counter):  # Step through the molecule list
                for turn in rotation["angles"]:  # Step through the angles to be performed
                    pbar.update(1)

                    # Rotated atoms
                    rotated = []

                    # Center the molecule on the center atom
                    new_rotamer = center_on_atom(rotamers[count], rotation["center"])

                    # Print the title of the file being rotated
                    # print(f"{new_rotamer.name}__a{rotation['ancr']}-c{rotation['center']}-{turn}deg")

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
                    rotamers.append(new_rotamer)

    # Perform file saving
    if save_com_files:
        for molecule in tqdm(rotamers, desc="Saving com files"):
            write_job_to_com(molecule.atoms, title=molecule.name, output=com_output)

    if save_images:
        kw = {'save': True, 'show': False, 'output': image_output}
        saving = partial(plot_structure, **kw)

        # Using pools results in ~4x speed performance boost when saving images.
        # 360 rotations in 00:02:45 -> 00:00:42
        with Pool() as pool:
            list(pool.imap(saving, tqdm(rotamers, desc="Saving images")))


if __name__ == "__main__":
    main()
    input("Calculating and saving complete. Press enter to close. ")
