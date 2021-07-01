"""
@Author: Jackson K Elowitt
@Start Date: May 14, 2021
@Contact: jkelowitt@protonmail.com
@Site: github.com/jkelowitt/Rotapy
@Version: v2.1

The end goal of this script is to be able to take in a Gaussian09 .log file,
and allow the user to rotate specific elements of the contained molecule, and
generate Gaussian09 input files for optimization on a supercomputer. Ideally,
multiple rotations may be performed simultaneously, such as in a nested for loop.

Main Changes remaining:
    - Theres a lot of faffing about switching between atom numbers and the atom itself.
        There's got to be a better way to handle this.
    - Try looking up mathematical graph based code to find a better way to code Molecule

"""
import winsound as ws
from copy import deepcopy
from math import ceil
from time import sleep

import PySimpleGUI as sg

from classes import Atom
from functions import (bonded_atom_search, center_on_atom, check_bonds, rotate_point_around_vector, save_structure,
                       show_structure)
from parsing import make_molecule_from_file, parsing_dict, write_job_to_com

settings = None

# Directly modified variables
file_types = tuple(("Valid Types", "*." + key) for key in parsing_dict)
rotamer_count = 1
tasks = []

# Entry cell width
ew = 7

choice_buttons = []
aToolTip = " When rotating an alcoholic hydrogen, this would be the carbon. "
cToolTip = " When rotating an alcoholic hydrogen, this would be the oxygen. "
angleToolTip = " When rotating an alcoholic hydrogen, this would be the size of the rotation steps. "


def parse_tasks(tasks):
    """Parse the string based task menu, and convert it to a rotation queue"""
    new_tasks = []
    for item in tasks:
        items = item.split(", ")
        ancr_atom_num = int(float(items[0]))
        center_atom_num = int(float(items[1]))
        angle = float(items[2])
        new_tasks.append([ancr_atom_num, center_atom_num, angle])

    return new_tasks


def make_queue_from_tasks(tasks, file):
    """Convert the task queue to the rotation queue"""
    rotation_queue = []
    base_compound = None

    for task in tasks:
        ancr_atom_num = task[0]
        center_atom_num = task[1]
        angle = task[2]

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


def generate_rotamers(base_compound, rotation_queue, window):
    """Generate the rotamers and return the rotamers in Molecule form"""
    # Grab the progress bar and update it as we go
    window["p_text"]("Rotating")

    # This is the starting position of the compound. Free from any rotations
    rotamers = [deepcopy(base_compound)]

    # Sort the rotation queue to increase calculation efficiency
    rotation_queue.sort(key=lambda x: len(x["rotatees"]), reverse=True)

    # Step through the rotation queue
    for i, rotation in enumerate(rotation_queue):

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

        # Update Progress bar
        window["pbar"].update_bar(i, rotamer_count)

    return rotamers


def settings_window(settings):
    """Open the settings window, and allow the user to specify settings"""

    # Default settings are used if the user doesn't change anything, or resets to default
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

    # If the user hasn't changed the settings before, use the default
    if settings is None:
        settings = default.copy()

    input_width = 20
    input_height = 10

    # Structure of right column of the settings window. Contains input boxes.
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

    # Structure of the left column of the settings window. Contains uneditable text.
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

    # Additional settings in the form of checkboxes.
    file_settings = sg.Col([
        [sg.CB("Sequentially Name Files", k="seq", default=settings["seq"])],
    ])

    # Final format of the settings window
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

    # Settings window event loop
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == "save_settings":
            sleep(0.1)
            settings["charge"] = values["charge"]
            settings["mul"] = values["mul"]
            settings["job"] = values["job"]
            settings["theory"] = values["theory"]
            settings["basis"] = values["basis"]
            settings["cores"] = values["cores"]
            settings["memory"] = values["memory"]
            settings["linda"] = values["linda"]
            settings["seq"] = values["seq"]

            # If the user saves the settings, close the settings window.
            break

        elif event == "reset":
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
    """Returns the formatted main window object"""

    # Data entry padding
    dep = (0, 0)

    # Data entry layouts
    anchor_layout = sg.Col([
        [sg.T("Anchor", tooltip=aToolTip)], [sg.In(s=(ew, 1), k="anchor", tooltip=aToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    center_layout = sg.Col([
        [sg.T("Center", tooltip=cToolTip)], [sg.In(s=(ew, 1), k="center", tooltip=cToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    angle_layout = sg.Col([
        [sg.T("Angle", tooltip=angleToolTip)], [sg.In(s=(ew, 1), k="angle", tooltip=angleToolTip)]
    ], pad=dep, element_justification="center", vertical_alignment="top")

    # Left column of the window layout.
    # Contains the data entry cells, the total rotamer count, task queue, and progress bar
    left_col = sg.Col([
        [anchor_layout, center_layout, angle_layout],
        [sg.Listbox(values=tasks, key="rotations", auto_size_text=True, size=(160 // ew, 14), no_scrollbar=True)],
        [sg.Text(f"Total Rotamers: {rotamer_count}", key="rot_count", size=(22, 1), font="Arial 10")],
        [sg.T("Progress Bar", k="p_text"), sg.Prog(max_value=1, size=(10, 10), k="pbar")],
    ], element_justification="center", vertical_alignment="top")

    # Button size
    bsz = (17, 1)
    bpad = (1, 1)

    # Right column of the main window
    # Contains all the file browser buttons, settings pop-up button, and execution button
    right_col = sg.Col([
        [sg.T("Import Molecule", )],
        [sg.I(k="input_file", s=(10, 1)),
         sg.FileBrowse(target="input_file", file_types=file_types, k="input_browse")],
        [sg.B(button_text="Show Molecule", k="show_plot", s=bsz, pad=bpad)],
        [sg.HSep()],
        [sg.T("Rotation Queue")],
        [sg.B(button_text='Add', key="add_save", s=bsz, pad=bpad)],
        [sg.B('Remove', s=bsz, pad=bpad)],
        [sg.HSep()],
        [sg.T("Output Settings")],
        [sg.T("Com Output"),
         sg.I(k="com_dir", s=(10, 1), tooltip=" If you don't want to save the com files, leave this blank. "),
         sg.FolderBrowse(target="com_dir", k="com_browse")],
        [sg.T("Img Output"),
         sg.I(k="img_dir", s=(10, 1), tooltip=" If you don't want to save the images, leave this blank. "),
         sg.FolderBrowse(target="img_dir", k="img_browse")],
        [sg.B("Change Output Settings", k="change_settings")],
        [sg.HSep()],
        [sg.B("Perform Calculations", k="execute")],

    ], element_justification="center", vertical_alignment="top")

    # Final layout of the main window
    layout = [
        [sg.Titlebar('Rotapy')],
        [left_col, right_col],
    ]

    return sg.Window('Rotapy', layout, keep_on_top=True, finalize=True)


window = make_main_window()

# Initialize rotamer count
for thing in parse_tasks(tasks):
    rotamer_count *= 360 // thing[-1]
window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

while True:
    # Wait till the user does something, then go into the condition loop
    event, values = window.Read()

    # Main event loop
    if event == "add_save":  # Add an item to the rotation queue
        try:
            a = values['anchor']
            c = values['center']
            d = values["angle"]

            if a == "" or c == "" or d == "":
                ws.MessageBeep()
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
            ws.MessageBeep()
            sg.popup_error("Entries into the rotation queue must be numerical", title="Rotation Queue Error",
                           keep_on_top=True)
            continue

        except AssertionError:
            ws.MessageBeep()
            sg.popup_error("The anchor atom and the center atom must not be the same atom",
                           title="Anchor Center Error", keep_on_top=True)
            continue

        if d > 360 or d < 0:
            ws.MessageBeep()
            sg.popup_error("The angle must be >0° and <360°", title="Angle Error", keep_on_top=True)
            continue

        tasks.append(f"{a}, {c}, {d}")

        # Calculate the number of rotations. Will always be an integer, hence floor_divide here
        rotamer_count *= 360 // d

        # Update rotation list
        window['rotations'].update(values=tasks)

        # Update rotation count
        window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

        # Clear inputs
        window['anchor'](value="")
        window['center'](value="")
        window['angle'](value="")

    elif event == "Remove":
        # Remove an item from the rotation queue
        try:
            v = values["rotations"][0].split(", ")
        except IndexError:
            ws.MessageBeep()
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

        window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

    elif event == "show_plot":
        """Reads the current import file, and shows the structure"""
        if not values["input_file"]:
            ws.MessageBeep()
            sg.popup_error("Cannot show plot until input file is entered.", title="Plotting Error", keep_on_top=True)
            continue

        show_plot(values)

    elif event == "execute":
        """Performs all the rotations in the rotation queue"""

        # Check that all the required values are present
        if not (values["input_file"] and tasks):
            ws.MessageBeep()
            sg.popup_error("Cannot perform calculations until both an item has been entered "
                           "into the rotation queue, and an input molecule has been selected.",
                           title="Exectution Error (1)", keep_on_top=True)
            continue

        elif not values["input_file"]:
            ws.MessageBeep()
            sg.popup_error("Cannot perform calculations until input file is entered.", title="Execution Error (2)",
                           keep_on_top=True)
            continue

        elif not tasks:
            ws.MessageBeep()
            sg.popup_error("Cannot perform calculations at least one task is entered.", title="Execution Error (3)",
                           keep_on_top=True)
            continue

        # If all goes well, perform the calculations
        file = values["input_file"]
        formatted_tasks = parse_tasks(tasks)
        rotation_queue, base_compound = make_queue_from_tasks(formatted_tasks, file)
        rotamers = generate_rotamers(base_compound, rotation_queue, window)

        # If the user doesn't change the settings, use these instead.
        # TODO this is a duplicate of the default in the settings_window function.
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
            window["p_text"]("Renaming Files")
            for i, rotamer in enumerate(rotamers):
                rotamer.name = f"{base_compound.name}_{i + 1}"
                window["pbar"].update_bar(i, rotamer_count)

        # Check for errors in the bonding
        errored_out = 0
        window["p_text"]("Checking ERR")
        for i, rotamer in enumerate(rotamers):
            if count := check_bonds(base_compound, rotamer):
                rotamer.name += f"_{count}ERR"
                errored_out += 1
            window["pbar"].update_bar(i, rotamer_count)

        # Alert the user to the presence of bad rotamers
        e_msg = f"{errored_out} rotamers had collisions while rotating. " \
                f"The files were marked with '_#ERR', where ## is " \
                f"the number of bonds different from the normal compound."

        if errored_out:
            # non_blocking is used so that this error message doesn't
            # prevent calculations from continuing when being run unattended.
            sg.popup_error(e_msg, title="Collisions Detected", keep_on_top=True, non_blocking=True)

        if output := values["com_dir"]:
            window["p_text"]("Writing COM")
            for i, molecule in enumerate(rotamers):
                write_job_to_com(molecule, title=molecule.name, directory=output, **settings)
                window["pbar"].update_bar(i)

        if output := values["img_dir"]:
            window["p_text"]("Writing IMG")
            for i, molecule in enumerate(rotamers):
                save_structure(molecule, directory=output)
                window["pbar"].update_bar(i)

        # Reset progress bar
        window["p_text"]("Progress Bar")
        window["pbar"].update_bar(0, rotamer_count)

        # Alert the user to the completion of the calculation
        ws.MessageBeep()
        sg.popup("Calculations Complete", keep_on_top=True)

    elif event == "change_settings":
        if x := settings_window(settings):
            settings = x.copy()

    elif event is None or event in ("Cancel", sg.WIN_CLOSED):
        break

window.Close()
