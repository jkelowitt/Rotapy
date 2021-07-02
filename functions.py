"""
Functions used to facilitate Rotapy.

plot_structure: Plots the structure of a given Molecule in 3d.

center_on_atom: Translates the coordinates of all atoms in the mo so that
                the selected atom is at the origin.

rotate_point_around_vector: Rotate a point around a selected vector by some number of degrees.
                            Uses a quaternion-ic rotation to avoid gimbal lock, and for ease of coding.

bonded_atom_search: Takes in a network and returns all nodes with a path to the starting atom,
                    excluding any blocked nodes.

verified_input: Verify that the user has input a value which can be converted to a specified type.
                This function will not return until the user has input something which can be converted
                to the type specified by 'verify'

"""

from os import getcwd, makedirs, path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyquaternion as pq

from classes import Atom, Molecule


def generate_figure(mo: Molecule):
    """
    Generate a figure of the molecule and bonds

    Parameters
    ----------
    mo: The molecule to be modeled

    Return
    ------
    fig: The matplotlib figure of the molecule.

    ax: The matplotlib axes of the molecule.
    """
    # Update the bond graph
    mo.make_bond_graph()

    # Make an empty, square, pseudo-3d figure
    dpi = 300
    pixels = 800
    fig = plt.figure(figsize=(pixels / dpi, pixels / dpi), dpi=dpi, tight_layout=True)
    ax = fig.add_subplot(projection="3d")

    # Hide axis planes and axis lines
    ax.set_axis_off()

    # Relative size scalar for the atoms.
    size = pixels / dpi * 50

    # Bonds are slightly see through and gray.
    alpha = 0.90
    line_color = '0.5'

    # Draw atoms as circles
    for num, a in enumerate(mo.atoms):
        x = a.pos[0]
        y = a.pos[1]
        z = a.pos[2]

        # Size scales with cov_radius
        ax.scatter(x, y, z, color=a.color, edgecolor=line_color, s=size * a.cov_radius)

        # Number the as.
        text_color = tuple(
            1 - a for a in a.color
        )  # The text color is the negative of the color of the a
        ax.text(x, y, z, num,
                zorder=100,
                color=text_color,
                ha="center",  # Horizontal Alignment
                va="center",  # Vertical Alignment
                fontsize=size * a.cov_radius / 16,
                )

    # Draw bonds as lines
    # This will draw duplicate lines on top of each other.
    # For example, if two atoms are bonded to each other,
    # there will be one bond from A to B, and another from B to A.
    # This may be a problem with very large molecules.
    for a in mo.bonds:
        x1 = a.pos[0]
        y1 = a.pos[1]
        z1 = a.pos[2]
        for bonded_a in mo.bonds[a]:
            x2 = bonded_a.pos[0]
            y2 = bonded_a.pos[1]
            z2 = bonded_a.pos[2]
            ax.plot((x1, x2), (y1, y2), (z1, z2), color=line_color, alpha=alpha)

    # One liner used to force the axes to be equal in step width.
    # For example, (0,1) may appear 20px in length, while (1,0) may appear 10px in length.
    # With this line, they will both be length 20.
    # I don't know how this works. I found it here:
    # https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
    ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")])

    return fig, ax


def show_structure(mo, title: str = None):
    """Show the structure of a molecule in an interactive plot"""
    matplotlib.use('tkAgg')

    fig, ax = generate_figure(mo)

    plt.title(title if title else mo.name, fontsize=5)
    plt.show()
    plt.close('all')


def save_structure(mo, title: str = None, directory: str = ""):
    """Save the structure of a molecule to a png file"""
    matplotlib.use('agg')

    fig, ax = generate_figure(mo)
    directory = f"{directory}/{title if title else mo.name}.png"

    plt.title(title if title else mo.name, fontsize=5)
    plt.savefig(directory)
    plt.close('all')


def center_on_atom(mo: Molecule, atom_number: int) -> Molecule:
    """
    Translates the coordinates of all atoms in the mo so that
        the selected atom is at the origin.

    Parameters
    ----------
    mo: The molecule to be shifted

    atom_number: The index of the atom in mo.atoms to be designated as the new center atom.

    Returns
    -------
    mo: A copy of the molecule shifted so that the atom indicated by atom_number is at (0,0,0).
    """
    selected_atom = mo.atoms[atom_number]
    dx = selected_atom.pos[0]
    dy = selected_atom.pos[1]
    dz = selected_atom.pos[2]

    new_atoms = list()
    for at in mo.atoms:
        new_position = (at.pos[0] - dx, at.pos[1] - dy, at.pos[2] - dz)
        new_at = Atom(name=at.name, pos=new_position)
        new_atoms.append(new_at)

    new_mo = Molecule(name=mo.name, atoms=new_atoms)

    return new_mo


def rotate_point_around_vector(point: tuple[float, float, float],
                               vector: tuple[float, float, float],
                               deg: float) -> tuple[float, float, float]:
    """
    Rotate a point around a selected vector by some number of degrees.
    Uses a quaternion-ic rotation to avoid gimbal lock, and for ease of coding.

    Parameters
    ----------
    point: The (x, y, z) point to be rotated.

    vector: The (x, y, z) vector to be used as a rotation axis. Center must be on the origin.
            Does not need to be unit size.

    deg: The angle to which the point will be rotated. Has units of degrees.

    Returns
    -------
    new_vec: A tuple containing the (x, y, z) cartesian coordinates of the rotated point
    """

    # This quaternion object will perform the rotation encoded into its initialization.
    # The rotation vector is of unit length.
    rotation_quaternion = pq.Quaternion(axis=vector, degrees=deg)

    # Perform the rotation
    new_vec = rotation_quaternion.rotate(point)
    return new_vec


def bonded_atom_search(molecule: Molecule, start: Atom, wall: list) -> list:
    """
    Takes in a network and returns all nodes with a path to the starting atom,
    excluding any blocked nodes.

    Parameters
    ----------
    molecule: A Molecule object

    start: An atom in the molecule that we want to find all the bonded atoms of

    wall:  A list of atoms which are not to be included in the search.

    Returns
    -------
    important: The list of atoms which have a path back to the start atom,
               excluding any paths through the wall atoms.
    """
    molecule.make_bond_graph()
    bonded = lambda x: molecule.bonds[x]
    important = [start]

    for atom in important:
        bonds = bonded(atom)
        for bond in bonds:
            if bond not in wall and (bond not in important):
                important.append(bond)

    return important


def verified_input(prompt: str = "", verify: type = int):
    """
    Verify that the user has input a value which can be converted to a specified type.
    This function will not return until the user has input something which can be converted
        to the type specified by 'verify'

    Parameters
    ----------
    prompt: The prompt for input for the user.
    verify: The type for the input to be returned as

    Returns
    -------
    data: Input from the user with the type, verify.

    """
    while True:
        data = input(prompt)

        try:
            data = verify(data)
            break
        except ValueError:
            print(f"Error: Must be of type {verify.__name__}")

    return data


def check_bonds(m1, m2):
    """Checks whether the two molecules have the same bonding structure
    If the bond count lists are not exactly equal, a collision occurred
    Only a very complex scenario would really defeat this detection method.
    A complex problem == complex solution, thus, I procrastinate.
    """
    # Update bond graphs
    m1.make_bond_graph()
    m2.make_bond_graph()

    # Makes a list of the bond counts for each atom
    m1struct = [len(m1.bonds[b]) for b in m1.bonds]
    m2struct = [len(m2.bonds[b]) for b in m2.bonds]

    m1struct.sort()
    m2struct.sort()

    bond_diff = abs(sum(m1struct) - sum(m2struct)) // 2
    return bond_diff


def randomly_orient(mo):
    """Randomly orient a given molecule in 3d space using uniform distributions"""
    mo = center_on_atom(mo, 0)
    new_mo = Molecule(name=mo.name, atoms=list())

    center_on_atom(mo, 0)

    rand_axis = tuple(np.random.uniform(size=(3,)))
    rand_deg = np.random.uniform() * 360

    for atom in mo.atoms:
        new_atom_pos = rotate_point_around_vector(point=atom.pos, vector=rand_axis, deg=rand_deg)
        new_atom = Atom(atom.name, new_atom_pos)
        new_mo.add_atom(new_atom)

    return new_mo


def make_output_folder(sub: str = "") -> str:
    """
    Makes a directory in the script location to output the downloaded files

    Parameters
    ----------
    sub: The name of the directory to be made.

    Returns
    -------
    dir_path: The directory pointing to :sub:

    """
    # Finds the current directory
    dir_path = getcwd()

    # Makes the path for the new folder
    dir_path = dir_path + fr"\{sub}"

    # If the folder doesn't exist, make it.
    if not path.exists(dir_path):
        try:
            makedirs(dir_path)
        except FileExistsError:
            # Sometimes this error pops when using threading or multiprocessing.
            pass
    return dir_path
