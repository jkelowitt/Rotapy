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
from os import path, makedirs

import matplotlib.pyplot as plt
import numpy as np
import pyquaternion as pq
from tqdm import tqdm

"""
Dictionaries containing useful data for elements.

Dictionaries:
    - cov_rads: Covalent Radii
    - at_masses: Atomic Masses

Taken from TMPChem's Computational Chemistry github page:
https://github.com/tmpchem/computational_chemistry
"""

# Margin of error on whether or not to make a bond
# Chosen arbitrarily.
bond_threshold = 1.2

# Covalent (or ionic) radii by atomic element (Angstroms) from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
cov_rads = {'H': 0.37, 'C': 0.77, 'O': 0.73, 'N': 0.75, 'F': 0.71,
            'P': 1.10, 'S': 1.03, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33, 'He': 0.30,
            'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B': 0.88, 'Na': 1.02,
            'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K': 1.38, 'Ca': 1.00, 'Sc': 0.75,
            'Ti': 0.86, 'V': 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
            'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
            'Se': 1.17, 'Kr': 1.03, 'X': 0.00}

# Relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {'H': 1.00794, 'C': 12.0107, 'O': 15.9994, 'N': 14.0067,
             'F': 18.9984, 'P': 30.9738, 'S': 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
             'I': 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
             'Be': 9.01218, 'B': 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
             'Si': 28.0855, 'K': 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
             'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
             'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
             'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X': 0.00000}

"""
Functions whose main purpose is to parse data from files or the user. 
"""


def yes_no(prompt: str) -> bool:
    """Returns True if Yes, False if No."""
    yes = ["y", "Y", "Yes", "1"]
    no = ["N", "n", "No", "0"]

    while True:
        done = input(f"{prompt} (y/n): ")
        if done in yes + no:
            break

    # If no
    if done in no:
        return False

    # If yes
    return True


def make_choice_list(choices: list[str], prompt: str = "Select one of the following (ex. 2):", ret_num: bool = False):
    """
    Prints a prompt and a list of choices for the user to select from.

    Parameters
    ----------

    choices: List of printable objects for the user to choose from.
    prompt: The question to ask the user prior to showing the list of choices.
    ret_num: Whether to return the index of the choice, rather than the actual item in choices.

    Returns
    -------
    If ret_num == False
        The item in the list which the person selected.
    If ret_num == True
        The index of the user's choice within the list.
    """
    print(prompt)
    for n, item in enumerate(choices):
        print(f"{n + 1}) {item}")

    # Ensure that the input is valid.
    while True:
        try:
            selection = int(input("Selection: "))

            # Check that the selection is within the range of the list.
            assert 1 <= selection <= len(choices)
            chosen = choices[selection - 1]
            break

        except ValueError:
            print("The selection must be an integer.")
        except AssertionError:
            print("The selection is not an option in the list.")
    if ret_num:
        return selection - 1

    return chosen


def make_choice_dict(choices: dict, prompt: str = "Select one of the following (ex. 2):"):
    """
    Prints a prompt and a list of choices for the user to select from.

    Parameters
    ----------

    choices: Dictionary of printable keys and values for the user to choose from.
    prompt: The question to ask the user prior to showing the list of choices.

    Returns
    -------
    If ret_num == False
        The item in the list which the person selected.
    If ret_num == True
        The index of the user's choice within the list.
    """
    print(prompt)
    for n, item in enumerate(choices):
        print(f"{item} = {choices[item]}")

    print("Exit")

    # Ensure that the input is valid.
    while True:
        selection = input("Selection: ")
        try:
            key = selection
            if selection in ("exit", "Exit"):
                return key

            value = choices[selection]  # Will error out if selection is not in choices.
            break
        except KeyError:
            print("That is not a valid selection. Please type out the full name of the item to be changed.")

    return key


def make_output_folder(sub: str = "") -> str:
    """Makes a directory in the script location to output the downloaded files"""
    # Finds the current directory
    dir_path = path.dirname(path.realpath(__file__))

    # Makes the path for the new folder
    dir_path = dir_path + fr"\{sub}"

    # If the folder doesn't exist, make it.
    if not path.exists(dir_path):
        makedirs(dir_path)
    return dir_path


def parse_opt_geom_from_log(file: str) -> list:
    """
    Given a .log file which contains an optimized geometry, extract the (x,y,z) cartesian coordinates.

    Returns the atoms in the format:
    [["Atom 1 name", X_coord, Y_coord, Z_coord],
    ["Atom 2 name", X_coord, Y_coord, Z_coord]
    ...
    ]
    """

    # Read the data from the file
    with open(file, "r+") as f:
        lines = f.readlines()  # Caution, files may be very /very/ large.

    # The cartesian data is the only data in the file which contains a \
    result_data = [line for line in lines if "\\" in line]

    # Combine the lines into a single line
    result_string = ""
    for line in result_data:
        result_string += line.replace("\n", "").replace(" ", "")

    # Split the data into the \'ed chunks, and remove everything which isn't the cartesian coordinates
    data = result_string.split("\\")[16:-11]

    molecule = []
    for entry in data:
        a = entry.split(",")

        name = a[0]
        x = float(a[1])
        y = float(a[2])
        z = float(a[3])

        new_entry = [name, x, y, z]
        molecule.append(new_entry)

    return molecule


def write_job_to_com(atoms: list,
                     title: str = "molecule_name",
                     charge: int = 0,
                     multiplicity: float = 1,
                     job: str = "Opt Freq",
                     theory: str = "B3LPY",
                     basis_set: str = "6-311G(2df,2p)",
                     cores: int = 8,
                     memory: str = "20gb",
                     linda: int = 1,
                     output: str = "") -> None:
    """
    Takes in a list of atoms and their cartesian coordinates such as in parse_opt_geom_from_log,
    and saves the coordinates to a .com file.

    The user may specify:
        The jobs to be performed
        The level of theory to use
        The basis set to run
        The number of cores to use
        The amount of memory to use
        Whether or not to use a checkpoint
        How many linda cores to use (set to 1 even if not being used)
        The output directory for the file
    """

    d = f"""\
%NProcShared={cores}
%NProcLinda={linda}
%mem={memory}
%Chk={title}.chk
#n {theory}/{basis_set} {job}

{title}

{charge} {multiplicity}
"""
    for a in atoms:
        name = str(a.name)
        x = str(a.pos[0])[:14].rjust(15)
        y = str(a.pos[1])[:14].rjust(15)
        z = str(a.pos[2])[:14].rjust(15)
        d += f"\n{name} {x} {y} {z}"

    directory = make_output_folder(output)

    with open(fr"{directory}\{title}.com", "w+") as file:
        file.write(d)


"""
Classes used for structuring out the data
"""


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


def center_on_atom(mo: Molecule,
                   atom_number: int) -> Molecule:
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


def distance(pt1: tuple[float, float, float],
             pt2: tuple[float, float, float]) -> float:
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

    # Get all the log files in the current directory and all subdirectories.
    files = glob("*.log")

    # Check that there are log files to be found.
    if not files:
        print("No log files found in the current directory or lower.")
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

        print(
            f"\nThe rotation of {len(rotate_nums)} atoms,\n"
            f"\tabout the {center_atom_num},{ancr_atom_num} atom bond axis,\n"
            f"\twith a total of {len(degrees) + 1}, {angle} degree rotations,\n"
            f"\thas been added to the queue of rotations.\n"
        )

        if not yes_no("Add more rotations"):
            break

    rotamers = [base_compound]

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

        settings = {"charge": "0",
                    "mul": "1",
                    "job": "Opt Freq",
                    "theory": "B3LPY",
                    "basis set": "6-311G(2df,2p)",
                    "cores": "8",
                    "mem": "20gb",
                    "linda": "1"}

        # Display default settings
        print("\nDefault Settings:")
        for item in settings:
            print(f"\t{item} = {settings[item]}")

        non_default = yes_no("Use the default settings")

        # Change default options if desired
        if not non_default:
            settings["Exit"] = ""
            while True:
                option = make_choice_dict(settings, prompt="\nWhich would you like to change (ex. job): ")

                if option in ("exit", "Exit"):
                    break

                value = input("What would you like to change it to (ex. opt): ")
                settings[option] = value
                print(f"Changed {option} to {value}\n")

    if save_images:
        image_output = input("What would you like to name the output directory for the image files: ")

    # Pause for preparation and alert the user to the download file count.
    _ = input(f"Press any key to save the {save_com_files * len(rotamers) + save_images * len(rotamers)} files.")

    # Perform file saving
    if save_com_files:
        for molecule in tqdm(rotamers, desc="Saving rotamer com files"):
            write_job_to_com(molecule.atoms, title=molecule.name, output=com_output)

    if save_images:
        for molecule in tqdm(rotamers, desc="Saving rotamer images"):
            molecule.plot_structure(save=True, show=False, output=image_output)


if __name__ == "__main__":
    main()
