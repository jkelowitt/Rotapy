"""
Classes used to facilitate Rotapy.

Atom: Data class containing all the pertinent information about an atom.
Molecule: A molecule class used to hold atoms and the bonds between those atoms.
"""

from dataclasses import dataclass, field
from math import sqrt

from numba import njit
from numpy import array, random

from data_dicts import atom_color, bond_threshold, cov_rads


# noinspection PyUnresolvedReferences
@dataclass(eq=True, unsafe_hash=True)
class Atom:
    """
    Data class containing all the pertinent information about an atom.
    eq = True: allows the comparison of atoms to each other
    unsafe_hash = True: allows the Atoms to be used as the keys of dictionaries

    If any parameter of the atom is changed, it will no longer work as the key of the dictionary.
        hence, 'unsafe' hash.

    Parameters
    ----------
    name: Atomic symbol of the atom

    pos: (x, y, z) position of the atom in global space

    color: Color of the atom. Defined in __post_init__ using a color dictionary

    cov_radius: Covalent Radius of the atom. Defines the distance at which the atom can bond.

    """

    name: str
    pos: tuple[float, float, float] = field(default=(0, 0, 0))

    def __post_init__(self):
        self.cov_radius: float = cov_rads[self.name]
        self.color: tuple[float, float, float] = atom_color.get(self.name, (1, 0, 1))  # Default color magenta


# noinspection PyUnresolvedReferences
class Molecule:
    """
    A molecule class used to hold atoms and the bonds between those atoms.

    Parameters
    ----------
    name: The name of the molecule.

    atoms: A list of atoms containing Atom objects which form the molecule

    bonds: A dictionary with Atoms as keys, and other Atoms as values.
           The key represents the current atom.
           The value represents every other atom which it is currently bonded to.

    Methods
    -------
    add_atom: Adds an atom to the molecule

    remove_atom: Removes an atom from the molecule

    replace_atom: Replaces an atom in the molecule, so that the new atom retains the index in
                  Molecule.atoms that the removed atom previously had.

    make_bond_graph: Generates self.bonds. Used after all atom adjustments.

    """

    def __init__(self, name: str, atoms: list):
        self.name = name
        self.atoms = atoms.copy()  # No mutability please
        self.bonds: dict = dict()

    def __repr__(self):
        return self.name

    def add_atom(self, other: Atom) -> None:
        """Add an atom to the molecule"""
        if isinstance(other, type(Atom)):
            raise TypeError(f"Must add an Atom to the molecule, not {type(other)}")

        self.atoms.append(other)

    def add_molecule(self, other, global_pos: tuple[float, float, float] = (0, 0, 0)):
        """Adds a group of atoms to the molecule"""
        a = center_on_atom(other, 0)
        for atom in a.atoms:
            new_pos = tuple(array(atom.pos) + array(global_pos))
            atom.pos = new_pos
            self.atoms.append(atom)

    def remove_atom(self, a) -> None:
        """Remove an atom from the molecule"""
        if isinstance(a, type(Atom)):
            raise TypeError(f"Must add an Atom to the molecule, not {type(a)}")

        self.atoms.remove(a)

    def replace_atom(self, old_num, new_atom) -> None:
        """
        Removes one atom and adds another in its place.
        This ensures that the numbering for the atoms in the molecule remain constant.
        Atom 2 will remain atom 2 even if atom 1 is replaced.
        """
        self.atoms[old_num] = new_atom

    def make_bond_graph(self) -> None:
        """
        Creates a dictionary with an entry for every atom in the molecule. For each atom, find
        the distance to every other atom. If the distance is within the sum of the covalent
        radii for each atom (within some margin of error), they are bonded. This bond is recorded
        in the dictionary entry for that atom.

        This is one of the slowest functions in this entire program, second only to the save_figure
        function. Any improvements to this function or its components yield massive returns.
        """

        self.bonds = dict()
        for a in self.atoms:
            new_bonds = []
            for other in self.atoms:
                max_bond_length = bond_threshold * (other.cov_radius + a.cov_radius)
                d = distance(other.pos, a.pos)
                if other != a and d <= max_bond_length:
                    new_bonds.append(other)

            self.bonds[a] = list(new_bonds).copy()


# noinspection PyUnresolvedReferences
@dataclass
class Sphere:
    """
    A dataclass for handling spheres

    Parameters
    ----------

    pos: A 3d tuple representing the 3d-cartesian coordinates of the center of the sphere
    radius: A float representing the radius of the sphere


    Functions
    ---------

    within: Returns True if a given point is within the sphere, False otherwise.

    """
    pos: tuple[float, float, float]
    radius: float

    def __repr__(self):
        return f"A sphere with radius {self.radius}, centered on {self.pos}"

    def within(self, point) -> bool:
        """
        Returns True if the given point is within the sphere

        Example
        -------
        >>> from classes import Sphere
        >>> s = Sphere(pos=(1, 1, 1), radius=5)
        >>> p = (0.5, 3, 2)
        >>> s.within(p)
        True

        Example
        -------

        >>> from classes import Sphere
        >>> s = Sphere(pos=(1, 1, 1), radius=5)
        >>> p = (200, 200, 200)
        >>> s.within(p)
        False

        """
        # Sphere position
        sx = self.pos[0]
        sy = self.pos[1]
        sz = self.pos[2]

        # Point position
        px = point[0]
        py = point[1]
        pz = point[2]

        return ((sx - px) ** 2) + ((sy - py) ** 2) + ((sz - pz) ** 2) <= (self.radius ** 2)

    def sample(self, source: str = "volume"):
        """
        Get a list of points normally distributed on the surface of the sphere.
        Has precision up to 15 sigfigs.

        Parameters
        ----------
        source: The source of the sampling. If == surface, the returned points
            will be on the surface of the sphere. If == anything else, the points
            returned will only be within the volume of the sphere.

        """
        point = (0, 0, 0)
        valid = False
        while not valid:
            # Get a random points in 3d space
            point = random.normal(0, 1, size=(3,))

            if source == "surface":
                # Normalize the point so that it lies on the unit sphere
                point /= sqrt(np.sum(square(point), axis=1)).reshape(point_count, 1)

            # Scale the sphere to be the same size as the desired sphere
            point *= self.radius

            # Shift the point into the correct position
            point += self.pos

            point = list(point)

            valid = self.within(point)

        return point


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


@njit
def distance(pt1: tuple[float, float, float], pt2: tuple[float, float, float]) -> float:
    """Returns the distance between two three-dimensional points"""
    return sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2 + (pt1[2] - pt2[2]) ** 2)
