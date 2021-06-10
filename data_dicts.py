"""
Dictionaries containing useful data for elements.

Dictionaries:
    - atom_color: Color representations
    - cov_rads: Covalent Radii
    - at_masses: Atomic Masses

Taken from TMPChem's Computational Chemistry github page:
https://github.com/tmpchem/computational_chemistry
"""

atom_color = {
    "C": (0, 0, 0),
    "H": (1, 1, 1),
    "O": (1, 0, 0),
    "S": (1, 1, 0),
    "Cl": (0, 1, 0),
}

# Margin of error on whether or not to make a bond
# Chosen arbitrarily.
bond_threshold = 1.2

# covalent (or ionic) radii by atomic element (Angstroms) from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
cov_rads = {'H': 0.37, 'C': 0.77, 'O': 0.73, 'N': 0.75, 'F': 0.71,
            'P': 1.10, 'S': 1.03, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33, 'He': 0.30,
            'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B': 0.88, 'Na': 1.02,
            'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K': 1.38, 'Ca': 1.00, 'Sc': 0.75,
            'Ti': 0.86, 'V': 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
            'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
            'Se': 1.17, 'Kr': 1.03, 'X': 0.00}

# relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {'H': 1.00794, 'C': 12.0107, 'O': 15.9994, 'N': 14.0067,
             'F': 18.9984, 'P': 30.9738, 'S': 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
             'I': 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
             'Be': 9.01218, 'B': 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
             'Si': 28.0855, 'K': 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
             'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
             'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
             'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X': 0.00000}
