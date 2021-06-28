"""
This is a working file when transferring over to a gui based input format
"""
import winsound as ws

import PySimpleGUI as sg

from functions import show_structure
from parsing import make_molecule_from_file, parsing_dict

file_types = tuple(("Valid Types", "*." + key) for key in parsing_dict)
rotamer_count = 1
tasks = []

# Entry cell width
ew = 7

choice_buttons = []
aToolTip = "When rotating an alcoholic hydrogen, this would be the carbon."
cToolTip = "When rotating an alcoholic hydrogen, this would be the oxygen."
angleToolTip = "When rotating an alcoholic hydrogen, this would be the hydrogen."

# Data entry cells and titles
# Pulled out so that the titles can be centered on the cell
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
    [sg.Listbox(values=tasks, key="rotations", auto_size_text=True, size=(160 // ew, 10), no_scrollbar=True,
                select_mode="LISTBOX_SELECT_MODE_EXTENDED")],
    [sg.Text(f"Total Rotamers: {rotamer_count}", key="rot_count", size=(22, 1), font="Arial 10")],
    [sg.Prog(max_value=1, size=(17, 1))],
], element_justification="center", vertical_alignment="top")

# Button size
bsz = (8, 1)
bpad = (1, 1)

right_col = sg.Col([
    [sg.T("Import Molecule", )],
    [sg.I(k="input_file", s=(10, 1)), sg.FileBrowse(target="input_file", file_types=file_types)],
    [sg.B(button_text="Show Molecule", k="show_plot", s=(17, 1), pad=bpad, )],
    [sg.HSep()],
    [sg.T("Rotation Queue")],
    [sg.B(button_text='Add', key="add_save", s=bsz, pad=bpad), sg.B('Delete', s=bsz, pad=bpad)],
    [sg.HSep()],
    [sg.T("Save Locations")],
    [sg.T("Com"), sg.I(k="com_dir", s=(10, 1)), sg.FolderBrowse(target="com_dir")],
    [sg.T("Img  "), sg.I(k="img_dir", s=(10, 1)), sg.FolderBrowse(target="img_dir")],
    [sg.HSep()],
    [sg.B("Perform Calculations", k="execute")],

], element_justification="center", vertical_alignment="top")

layout = [
    [sg.Titlebar('Rotapy')],
    [left_col, right_col],
]


def show_plot():
    ...


def run_calculations():
    ...


window = sg.Window('Rotapy', layout)
while True:  # Event Loop
    event, values = window.Read()
    try:

        if event == "add_save":  # Add an item to the rotation queue
            a = int(values['anchor'])
            c = int(values['center'])
            d = float(values["angle"])

            tasks.append(f"{a}, {c}, {d}")

            # Calculate the number of rotations. Will always be an integer, hence floor_divide here
            rotamer_count *= 360 // d

            window['rotations'].update(values=tasks)

            rot_string = str(f"Total Rotamers: {rotamer_count}")
            window['rot_count']("Total Rotamers: {}".format(int(rotamer_count)))

        elif event == "Delete":  # Remove an item from the rotation queue
            v = values["rotations"][0].split(", ")

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
            """Reads the current import file, and shows the structure in a non-blocking way"""
            ws.Beep(1000, 100)
            show_plot()

        elif event == "execute":
            """Performs all the rotations in the rotation queue"""
            sg.popup('Some rotamers were found to have errors. They are marked with "ERR" in their name.',
                     title="Errors Found")
            run_calculations()

        elif event is None or event in ("Cancel", sg.WIN_CLOSED):
            break

    # Problems? Scream!
    except IndexError:
        ws.MessageBeep()

    except ValueError:
        ws.MessageBeep()

window.Close()
