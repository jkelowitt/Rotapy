import matplotlib
import matplotlib.pyplot as plt
import pyquaternion as pq

from classes import Atom, Molecule
from parsing import make_output_folder


def plot_structure(mo, title: str = None, save: bool = False, show: bool = True, output: str = "") -> None:
    """
    Plots the structure of the Molecule in 3d.
    The atoms are colored according to the color_dict in __post_init__.
    The atoms are sized proportional to the atom's cov_radius.
    The atoms are numbered by the order in which they are added to the molecule.
    Bonds are placed according to mo.bond_graph.
    Double and triple bonds are not shown.
    """

    # Prevent memory leak when saving many figures.
    if save:
        matplotlib.use('agg')
    else:
        matplotlib.use('tkAgg')

    # Make an empty 3d figure
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Hide axis planes and axis lines
    ax.set_axis_off()

    dpi = 200
    fig.set_dpi(dpi)
    fig.set_size_inches(3, 3)

    # Size = 300 @ 140 dpi
    size = dpi / 2

    # Draw atoms as circles
    for num, a in enumerate(mo.atoms):
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
                fontsize=size * a.cov_radius / 10,
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
            ax.plot((x1, x2), (y1, y2), (z1, z2), color=(0.5, 0.5, 0.5))

    # One liner used to force the axes to be equal in step width.
    # For example, (0,1) may appear 20px in length, while (1,0) may appear 10px in length.
    # With this line, they will both be length 20.
    # I don't know how this works. I found it here:
    # https://github.com/matplotlib/matplotlib/issues/17172#issuecomment-830139107
    ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")])

    if save:
        directory = f"{make_output_folder(output)}/{title if title else mo.name}.png"
        plt.savefig(directory)

    if show:
        plt.title(title if title else mo.name)
        plt.show()

    # Try to remove the figures from memory.
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
