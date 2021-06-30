"""
This is a working file when transferring over to a gui based input format
"""
import winsound as ws
from copy import deepcopy
from math import ceil

import PySimpleGUI as sg
from tqdm import tqdm

from classes import Atom
from functions import (bonded_atom_search, center_on_atom, check_bonds, rotate_point_around_vector, save_structure,
                       show_structure)
from parsing import make_molecule_from_file, parsing_dict, write_job_to_com

settings = None

# Directly modified variables
file_types = tuple(("Valid Types", "*." + key) for key in parsing_dict)
rotamer_count = 1
tasks = ["7, 17, 1"]

# Entry cell width
ew = 7

choice_buttons = []
aToolTip = " When rotating an alcoholic hydrogen, this would be the carbon. "
cToolTip = " When rotating an alcoholic hydrogen, this would be the oxygen. "
angleToolTip = " When rotating an alcoholic hydrogen, this would be the size of the rotation steps. "


def confirm_sound():
    """Play confirmation sound"""
    ws.Beep(1000, 100)


def deny_sound():
    """Play denial sound"""
    ws.Beep(250, 100)


def error_sound():
    """Play denial sound"""
    ws.MessageBeep()


def parse_tasks(queue, file):
    """Parse the string based task menu, and convert it to a rotation queue"""
    rotation_queue = []
    base_compound = None

    for item in queue:
        # Get user input numbers
        items = item.split(", ")
        ancr_atom_num = int(float(items[0]))
        center_atom_num = int(float(items[1]))
        angle = float(items[2])

        # Get all the angles from 0 degree rotation
        count = ceil((360 - angle) / angle)
        degrees = [round((n + 1) * angle, 5) for n in range(count)]

        # Obtain the physical atoms
        base_compound = make_molecule_from_file(file)
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

        rotation_queue.append({
            "center": center_atom_num,
            "ancr": ancr_atom_num,
            "rotatees": rotate_nums,
            "angles": degrees,
        })

    return rotation_queue, base_compound


def show_plot(v):
    """Show a plot of the input file"""
    molecule = make_molecule_from_file(v["input_file"])
    show_structure(molecule, title=molecule.name)


def generate_rotamers(base_compound, rotation_queue):
    """Generate the rotamers and return the rotamers in Molecule form"""
    rotamers = [deepcopy(base_compound)]

    # Sort the rotation queue to increase calculation efficiency
    rotation_queue.sort(key=lambda x: len(x["rotatees"]), reverse=True)

    # Step through the rotation queue
    for rotation in rotation_queue:

        # Counter holds the number of molecules to have rotamers made from.
        # This is done to prevent performing rotations twice on the same rotamer.
        counter = len(rotamers)

        # Step through the molecule list
        for count in range(counter):

            # Step through the angles to be performed
            for turn in rotation["angles"]:
                # Rotated atoms
                rotated = []

                # The molecule being rotated
                originator = rotamers[count]

                # Center the molecule on the center atom
                new_rotamer = center_on_atom(originator, rotation["center"])

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

                new_rotamer.name += f"_a{rotation['ancr']}c{rotation['center']}d{turn}"
                new_rotamer.make_bond_graph()
                rotamers.append(new_rotamer)
    return rotamers


def settings_window(settings):
    """Open the settings window"""
    default = {
        "charge": "0",
        "mul": "1",
        "job": "Opt Freq",
        "theory": "B3LYP",
        "basis": "6-311G(2df,2p)",
        "cores": "8",
        "memory": "20gb",
        "linda": "1",
        "seq": True
    }

    if settings is None:
        settings = default.copy()

    input_width = 20
    input_height = 10

    settings_right = sg.Col(
        [[sg.I(settings["charge"], k="charge", s=(input_width, input_height))],
         [sg.I(settings["mul"], k="mul", s=(input_width, input_height))],
         [sg.I(settings["job"], k="job", s=(input_width, input_height))],
         [sg.I(settings["theory"], k="theory", s=(input_width, input_height))],
         [sg.I(settings["basis"], k="basis", s=(input_width, input_height))],
         [sg.I(settings["cores"], k="cores", s=(input_width, input_height))],
         [sg.I(settings["memory"], k="memory", s=(input_width, input_height))],
         [sg.I(settings["linda"], k="linda", s=(input_width, input_height))],
         ])

    settings_left = sg.Col([
        [sg.T("charge", k="tch")],
        [sg.T("mul", k="tm")],
        [sg.T("job", k="tj")],
        [sg.T("theory", k="tt")],
        [sg.T("basis", k="tb")],
        [sg.T("cores", k="tc")],
        [sg.T("memory", k="tmem")],
        [sg.T("linda", k="tl")],
    ])

    file_settings = sg.Col([
        [sg.CB("Sequentially Name Files", k="seq", default=settings["seq"])],
        [],
        [],
        [],
        [],
    ])

    settings_layout = [
        [sg.Titlebar('Rotapy')],
        [sg.T("Job Settings:")],
        [settings_left, settings_right],
        [sg.HSep()],
        [sg.T("File Settings:")],
        [file_settings],
        [sg.HSep()],
        [sg.B("Save", k="save_settings"), sg.B("Reset to Default", k="reset")]
    ]
    window = sg.Window("Output Settings", settings_layout, keep_on_top=True)

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == "save_settings":
            confirm_sound()
            settings["charge"] = values["charge"]
            settings["mul"] = values["mul"]
            settings["job"] = values["job"]
            settings["theory"] = values["theory"]
            settings["basis"] = values["basis"]
            settings["cores"] = values["cores"]
            settings["memory"] = values["memory"]
            settings["linda"] = values["linda"]
            settings["seq"] = values["seq"]

            break

        elif event == "reset":
            deny_sound()
            window["charge"](default["charge"])
            window["mul"](default["mul"])
            window["job"](default["job"])
            window["theory"](default["theory"])
            window["basis"](default["basis"])
            window["cores"](default["cores"])
            window["memory"](default["memory"])
            window["linda"](default["linda"])
            window["seq"](default["seq"])

    window.close()
    return settings


def make_main_window():
    """Data entry cells and titles. Pulled out so that the titles can be centered on the cell"""
    dep = (0, 0)
    anchor_layout = sg.Col([
        [sg.T("Anchor", tooltip=aToolTip)], [sg.In(s=(ew, 1), k="anchor", tooltip=aToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    center_layout = sg.Col([
        [sg.T("Center", tooltip=cToolTip)], [sg.In(s=(ew, 1), k="center", tooltip=cToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    angle_layout = sg.Col([
        [sg.T("Angle", tooltip=angleToolTip)], [sg.In(s=(ew, 1), k="angle", tooltip=angleToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    left_col = sg.Col([
        [anchor_layout, center_layout, angle_layout],
        [sg.Listbox(values=tasks, key="rotations", auto_size_text=True, size=(160 // ew, 14), no_scrollbar=True)],
        [sg.Text(f"Total Rotamers: {rotamer_count}", key="rot_count", size=(22, 1), font="Arial 10")],
        [sg.Prog(max_value=1, size=(17, 1))],
    ], element_justification="center", vertical_alignment="top")

    # Button size
    bsz = (17, 1)
    bpad = (1, 1)

    right_col = sg.Col([
        [sg.T("Import Molecule", )],
        [sg.I(k="input_file", s=(10, 1), default_text=r"F:\Coding Projects\Python\Rotapy\example data\terpineol4.log"),
         sg.FileBrowse(target="input_file", file_types=file_types)],
        [sg.B(button_text="Show Molecule", k="show_plot", s=bsz, pad=bpad, )],
        [sg.HSep()],
        [sg.T("Rotation Queue")],
        [sg.B(button_text='Add', key="add_save", s=bsz, pad=bpad)],
        [sg.B('Remove', s=bsz, pad=bpad)],
        [sg.HSep()],
        [sg.T("Output Settings")],
        [sg.T("Com Output"),
         sg.I(k="com_dir", s=(10, 1), tooltip=" If you don't want to save the com files, leave this blank. ",
              default_text=r"F:\Coding Projects\Python\Rotapy\example data\com_output"),
         sg.FolderBrowse(target="com_dir")],
        [sg.T("Img Output"),
         sg.I(k="img_dir", s=(10, 1), tooltip=" If you don't want to save the images, leave this blank. ",
              default_text=r"F:\Coding Projects\Python\Rotapy\example data\img_output"),
         sg.FolderBrowse(target="img_dir")],
        [sg.B("Change Output Settings", k="change_settings")],
        [sg.HSep()],
        [sg.B("Perform Calculations", k="execute")],

    ], element_justification="center", vertical_alignment="top")

    layout = [
        [sg.Titlebar('Rotapy')],
        [left_col, right_col],
    ]

    return sg.Window('Rotapy', layout, keep_on_top=True)


window = make_main_window()
while True:
    """"""

    event, values = window.Read()
    if event == "add_save":  # Add an item to the rotation queue
        try:
            a = values['anchor']
            c = values['center']
            d = values["angle"]

            if a == "" or c == "" or d == "":
                error_sound()
                sg.popup_error(
                    "One or more of the inputs is empty. Please enter a selection to all"
                    " input cells before adding to the rotation queue.",
                    title="Empty Cell Error", keep_on_top=True)
                continue

            # int(float()) because 1.5 cannot be turned into an int while it is still a string
            a = int(float(a))
            c = int(float(c))
            d = float(d)

            assert a != c

        except ValueError:
            error_sound()
            sg.popup_error("Entries into the rotation queue must be numerical", title="Rotation Queue Error",
                           keep_on_top=True)
            continue

        except AssertionError:
            error_sound()
            sg.popup_error("The anchor atom and the center atom must not be the same atom",
                           title="Anchor Center Error", keep_on_top=True)
            continue

        if d > 360:
            error_sound()
            sg.popup_error("The angle must be less than 360°", title="Angle Error", keep_on_top=True)
            continue

        tasks.append(f"{a}, {c}, {d}")

        # Calculate the number of rotations. Will always be an integer, hence floor_divide here
        rotamer_count *= 360 // d

        # Update rotation list
        window['rotations'].update(values=tasks)

        # Update rotation count
        rot_string = str(f"Total Rotamers: {rotamer_count}")
        window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

        # Clear inputs
        window['anchor'](value="")
        window['center'](value="")
        window['angle'](value="")

    elif event == "Remove":  # Remove an item from the rotation queue
        try:
            v = values["rotations"][0].split(", ")
        except IndexError:
            error_sound()
            sg.popup_error("Click on one of the items in the rotation queue in order to remove it.",
                           title="Remove Error", keep_on_top=True)
            continue

        a = int(v[0])
        c = int(v[1])
        d = float(v[2])

        tasks_remove = f"{a}, {c}, {d}"
        tasks.remove(values["rotations"][0])

        # Calculate the number of rotations. Will always be an integer, hence floor_divide here
        rotamer_count /= 360 // d

        window['rotations'].update(values=tasks)

        rot_string = str(f"Total Rotamers: {rotamer_count}")
        window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

    elif event == "show_plot":
        """Reads the current import file, and shows the structure"""
        if not values["input_file"]:
            error_sound()
            sg.popup_error("Cannot show plot until input file is entered.", title="Plotting Error", keep_on_top=True)
            continue

        show_plot(values)

    elif event == "execute":
        """Performs all the rotations in the rotation queue"""

        # Check that all the required values are present
        if not (values["input_file"] and tasks):
            error_sound()
            sg.popup_error("Cannot perform calculations until both an item has been entered "
                           "into the rotation queue, and an input molecule has been selected.",
                           title="Exectution Error (1)", keep_on_top=True)
            continue

        elif not values["input_file"]:
            error_sound()
            sg.popup_error("Cannot perform calculations until input file is entered.", title="Execution Error (2)",
                           keep_on_top=True)
            continue

        elif not tasks:
            error_sound()
            sg.popup_error("Cannot perform calculations at least one task is entered.", title="Execution Error (3)",
                           keep_on_top=True)
            continue

        # If all goes well, perform the calculations
        file = values["input_file"]
        rotation_queue, base_compound = parse_tasks(tasks, file)
        rotamers = generate_rotamers(base_compound, rotation_queue)

        if settings is None:
            settings = {
                "charge": "0",
                "mul": "1",
                "job": "Opt Freq",
                "theory": "B3LYP",
                "basis": "6-311G(2df,2p)",
                "cores": "8",
                "memory": "20gb",
                "linda": "1",
                "seq": True
            }

        # Set the naming scheme to sequential if requested
        if settings["seq"]:
            for i, rotamer in enumerate(rotamers):
                rotamer.name = f"{base_compound.name}_{i + 1}"

        # Check for errors in the bonding
        errored_out = 0
        for rotamer in rotamers:
            if count := check_bonds(base_compound, rotamer):
                rotamer.name += f"_{count}ERR"
                errored_out += 1

        # Alert the user to the presence of bad rotamers
        e_msg = f"{errored_out} rotamers had collisions while rotating. " \
                f"The files were marked with '_#ERR', where ## is " \
                f"the number of bonds different from the normal compound."

        if errored_out:
            # non_blocking is used so that this error message doesn't
            # prevent calculations from continuing when being run unattended.
            sg.popup_error(e_msg, title="Collisions Detected", keep_on_top=True, non_blocking=True)

        if output := values["com_dir"]:
            for molecule in tqdm(rotamers, desc="Saving com files", dynamic_ncols=True):
                write_job_to_com(molecule, title=molecule.name, directory=output, **settings)

        if output := values["img_dir"]:
            for rotamer in tqdm(rotamers, desc="Saving img files", dynamic_ncols=True):
                save_structure(rotamer, directory=output)

    elif event == "change_settings":
        if x := settings_window(settings):
            settings = x.copy()

    elif event is None or event in ("Cancel", sg.WIN_CLOSED):
        break

window.Close()
