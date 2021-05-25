from os import path, makedirs


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
