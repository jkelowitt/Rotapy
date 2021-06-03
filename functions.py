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

import matplotlib
import matplotlib.pyplot as plt
import pyquaternion as pq

from classes import Atom, Molecule
from parsing import make_output_folder


def generate_figure(mo):
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

    # Make an empty 3d figure
    dpi = 300
    pixels = 800
    fig = plt.figure(figsize=(pixels / dpi, pixels / dpi), dpi=dpi, tight_layout=True)
    ax = fig.add_subplot(projection="3d")

    # Hide axis planes and axis lines
    ax.set_axis_off()

    # Relative size scalar for the atoms.
    size = pixels / dpi * 50

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


def save_structure(mo, title: str = None, output: str = ""):
    """Save the structure of a molecule to a png file"""
    matplotlib.use('agg')

    fig, ax = generate_figure(mo)
    directory = f"{make_output_folder(output)}/{title if title else mo.name}.png"

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

    new_mo = Molecule(name=mo.name, atoms=list().copy())
    for at in mo.atoms:
        new_position = (at.pos[0] - dx, at.pos[1] - dy, at.pos[2] - dz)
        new_at = Atom(name=at.name, pos=new_position)
        new_mo.add_atom(new_at)

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
    rotation_quaternion = pq.Quaternion(axis=vector, degrees=deg).unit

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
