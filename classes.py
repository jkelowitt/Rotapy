"""
Classes used to facilitate Rotapy.

Atom: Data class containing all the pertinent information about an atom.
Molecule: A molecule class used to hold atoms and the bonds between those atoms.
"""

from dataclasses import dataclass, field

from numpy import array
from numpy.linalg import norm

from data_dicts import cov_rads, bond_threshold


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
        self.bonds: dict = dict().copy()
        self.make_bond_graph()

    # @staticmethod
    # def distance(pt1: tuple[float, float, float], pt2: tuple[float, float, float]) -> float:
    #     """Returns the distance between two points (x, y, z) in 3D"""
    #     return sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2 + (pt1[2] - pt2[2]) ** 2)

    def add_atom(self, other: Atom) -> None:
        """Add an atom to the molecule"""
        if isinstance(other, type(Atom)):
            raise TypeError(f"Must add an Atom to the molecule, not {type(other)}")

        self.atoms.append(other)
        self.make_bond_graph()

    def remove_atom(self, a) -> None:
        """Remove an atom from the molecule"""
        self.atoms.remove(a)
        self.make_bond_graph()

    def replace_atom(self, old_num, new_atom) -> None:
        """
        Removes one atom and adds another in its place.
        This ensures that the numbering for the atoms in the molecule remain constant.
        Atom 2 will remain atom 2 even if atom 1 is replaced.
        """
        self.atoms[old_num] = new_atom
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
                distance = norm(array(other.pos) - array(a.pos))
                if other != a and distance <= max_bond_length:
                    new_bonds.append(other)

            self.bonds[a] = list(new_bonds).copy()