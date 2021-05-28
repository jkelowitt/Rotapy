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

import sys
from dataclasses import dataclass, field
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pyquaternion as pq
from tqdm import tqdm

from data_dicts import cov_rads, bond_threshold
from parsing import (
    make_output_folder,
    make_choice_list,
    yes_no,
    parse_opt_geom_from_log,
    make_choice_dict,
    write_job_to_com,
)


@dataclass(eq=True, unsafe_hash=True)
class Atom:
    """
    Data class containing all the pertinent information about an atom.
    """

    name: str
    pos: tuple[float, float, float] = field(default=(0, 0, 0))

    def __post_init__(self):
        colors = {
            "C": (0, 0, 0),
            "H": (1, 1, 1),
            "O": (1, 0, 0),
            "S": (1, 1, 0),
            "Cl": (0, 1, 0),
        }

        self.cov_radius: float = cov_rads[self.name]
        self.color: tuple[float, float, float] = colors[self.name]


class Molecule:
    """
    A molecule class used to hold atoms and the bonds between those atoms.
    """

    def __init__(self, name: str, atoms: list):
        self.name = name
        self.atoms = atoms.copy()  # No mutability please
        self.bonds: dict = dict().copy()
        self.make_bond_graph()

    def add_atom(self, other: Atom) -> None:
        """Add an atom to the molecule"""
        if isinstance(other, type(Atom)):
            raise TypeError(f"Must add an Atom to the molecule, not {type(other)}")

        self.atoms.append(other)
        self.make_bond_graph()

    def remove_atom(self, a):
        """Remove an atom from the molecule"""
        self.atoms.remove(a)
        self.make_bond_graph()

    def make_bond_graph(self) -> None:
        """
        Creates a dictionary with an entry for every atom in the molecule. For each atom, find
        the distance to every other atom. If the distance is within the sum of the covalent
        radii for each atom (within some margin of error), they are bonded. This bondis recorded
        in the dictionary entry for that atom.
        """

        self.bonds = dict()
        for a in self.atoms:
            new_bonds = []
            for other in self.atoms:
                max_bond_length = bond_threshold * (other.cov_radius + a.cov_radius)
                if other != a and distance(other.pos, a.pos) <= max_bond_length:
                    new_bonds.append(other)

            self.bonds[a] = list(new_bonds).copy()

    def plot_structure(self, title=None, save=False, show=True, output="") -> None:
        """
        Plots the structure of the Molecule in 3d.
        The atoms are colored according to the color_dict in __post_init__.
        The atoms are sized proportional to the atom's cov_radius.
        The atoms are numbered by the order in which they are added to the molecule.
        Bonds are placed according to self.bond_graph.
        Double and triple bonds are not shown.
        """

        # Regenerate the bonds in case any new atoms have been added.
        self.make_bond_graph()

        # Make an empty 3d figure
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        size = 300

        # Draw atoms as circles
        for num, a in enumerate(self.atoms):
            x = a.pos[0]
            y = a.pos[1]
            z = a.pos[2]

            # Size scales with cov_radius
            ax.scatter(x, y, z, color=a.color, edgecolor="k", s=size * a.cov_radius)

            # Number the as.
            text_color = tuple(
                1 - a for a in a.color
            )  # The text color is the negative of the color of the a
            ax.text(x, y, z, num,
                    zorder=100,
                    color=text_color,
                    ha="center",  # Horizontal Alignment
                    va="center",  # Vertical Alignment
                    fontsize=size * a.cov_radius / 20,
                    )

        # Draw bonds as lines
        # This will draw duplicate lines on top of each other.
        # For example, if two atoms are bonded to each other,
        # there will be one bond from A to B, and another from B to A.
        # This may be a problem with very large molecules.
        for a in self.bonds:
            x1 = a.pos[0]
            y1 = a.pos[1]
            z1 = a.pos[2]
            for bonded_a in self.bonds[a]:
                x2 = bonded_a.pos[0]
                y2 = bonded_a.pos[1]
                z2 = bonded_a.pos[2]
                ax.plot((x1, x2), (y1, y2), (z1, z2), color=(0.5, 0.5, 0.5))

        # One liner used to force the axes to be equal in step width.
        # For example, (0,1) may appear 20px in length, while (1,0) may appear 10px in length.
        # With this line, they will both be length 20.
        # I don't know how this works. I found it here:
        # https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
        ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")])

        if save:
            directory = make_output_folder(output)
            if title:
                plt.savefig(f"{directory}/{title}.png")
            else:
                plt.savefig(f"{directory}/{self.name}.png")
        if show:
            if title:
                plt.title(title)
            else:
                plt.title(self.name)
            plt.show()

        plt.close()

    def replace_atom(self, old_num, new_atom):
        """
        Removes one atom and adds another in its place.
        This ensures that the numbering for the atoms in the molecule remain constant.
        Atom 2 will remain atom 2 even if atom 1 is replaced.
        """
        self.atoms[old_num] = new_atom
        self.make_bond_graph()


def center_on_atom(mo: Molecule, atom_number: int) -> Molecule:
    """
    Translates the coordinates of all atoms in the mo so that
        the selected atom is at the origin.

    :return: A shifted copy of the molecule.
    """
    selected_atom = mo.atoms[atom_number]
    dx = selected_atom.pos[0]
    dy = selected_atom.pos[1]
    dz = selected_atom.pos[2]

    new_mo = Molecule(name=mo.name, atoms=list().copy())
    for at in mo.atoms:
        new_position = (at.pos[0] - dx, at.pos[1] - dy, at.pos[2] - dz)
        new_at = Atom(name=at.name, pos=new_position)
        new_mo.add_atom(new_at)

    return new_mo


def distance(pt1: tuple[float, float, float], pt2: tuple[float, float, float]) -> float:
    """Returns the distance between two points (x, y, z) in 3D"""
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2 + (pt1[2] - pt2[2]) ** 2)


def rotate_point_around_vector(point: tuple[float, float, float],
                               vector: tuple[float, float, float],
                               deg: float) -> tuple[float, float, float]:
    """
    Rotate a point around a selected vector by some number of degrees.
    Uses a quaternionic rotation to avoid gimbal lock, and for ease of coding.
    :param point: The point to be rotated.
    :param vector: The vector to be used as a rotation axis. Center must be on the origin.
    :param deg: The angle to which the point will be rotated. Has units of degrees.
    :return: A tuple containing the (x, y, z) cartesian coordinates of the rotated point
    """

    # This quaternion object will perform the rotation encoded into its initialization.
    # The rotation vector is of unit length.
    rotation_quaternion = pq.Quaternion(axis=vector, degrees=deg).unit

    # Perform the rotation
    new_vec = rotation_quaternion.rotate(point)
    return new_vec


def bonded_atom_search(molecule, start, wall):
    """
    Takes in a network and returns all nodes with a path to the starting atom,
    excluding any blocked nodes.
    """
    bonded = lambda x: molecule.bonds[x]
    important = [start]

    for atom in important:
        bonds = bonded(atom)
        for bond in bonds:
            if bond != wall and (bond not in important):
                important.append(bond)

    return important


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
        base_compound.plot_structure()

    while True:
        ancr_atom_num = int(input("Which atom is the anchor atom (ex. 7): "))
        center_atom_num = int(input("Which atom is the center atom (ex. 27): "))

        base_compound = center_on_atom(base_compound, center_atom_num)
        ancr_atom = base_compound.atoms[ancr_atom_num]
        center_atom = base_compound.atoms[center_atom_num]

        # All things attached to the center atom will be rotated.
        # The first atom is the center atom. Remove it to prevent numbering confusion.
        rotate_atoms = bonded_atom_search(base_compound, start=center_atom, wall=ancr_atom)[1:]

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
        degrees = np.arange(angle, 360, angle)

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

        non_default = yes_no("Use the default settings")

        # Change default options if desired
        if not non_default:
            settings["Exit"] = ""
            while True:
                option = make_choice_dict(
                    settings, prompt="\nWhich would you like to change (ex. job): "
                )

                if option in ("exit", "Exit"):
                    break

                value = input("What would you like to change it to (ex. opt): ")
                settings[option] = value
                print(f"Changed {option} to {value}\n")

    if save_images:
        image_output = input(
            "What would you like to name the output directory for the image files: "
        )

    # Pause for preparation and alert the user to the download file count.
    _ = input(
        f"Press any key to save the {save_com_files * len(rotamers) + save_images * len(rotamers)} files."
    )

    # Perform rotation calculations
    with tqdm(total=rotation_count, desc="Performing rotation calculations") as pbar:
        for rotation in rotation_queue:  # Step through the rotation queue
            counter = len(rotamers)  # Counter holds the number of molecules to have rotamers made from.
            for count in range(counter):  # Step through the molecule list
                for turn in rotation["angles"]:  # Step through the angles to be performed

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
                    pbar.update(1)

    # Perform file saving
    if save_com_files:
        for molecule in tqdm(rotamers, desc="Saving rotamer com files"):
            write_job_to_com(molecule.atoms, title=molecule.name, output=com_output)

    if save_images:
        for molecule in tqdm(rotamers, desc="Saving rotamer images"):
            molecule.plot_structure(save=True, show=False, output=image_output)


if __name__ == "__main__":
    main()
