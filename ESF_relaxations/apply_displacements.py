import sys

# Check if the user provided current_step as a command-line argument
if len(sys.argv) < 2:
    print("Error: Please provide the current step as a command-line argument.")
    sys.exit(1)

# Get current_step from command-line argument
current_step = int(sys.argv[1])

QE_files_to_modify = ["1-scf/scf.in",
                      "2-wfn/bands.in",
                      "3-wfnq/bands.in",
                      "4-wfn_co/bands.in"]

displacements_file = f"../step_{current_step-1}/11-excited_state_forces/displacements_Newton_method.dat"

def read_qe_input(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atomic_positions_start = None
    atomic_positions_end = None
    atomic_positions = []

    for i, line in enumerate(lines):
        if line.strip().startswith('ATOMIC_POSITIONS'):
            atomic_positions_start = i + 1
        elif atomic_positions_start and line.strip() == '':
            atomic_positions_end = i
            break
        elif atomic_positions_start:
            atomic_positions.append(line)

    return lines, atomic_positions, atomic_positions_start, atomic_positions_end

def read_displacements(file_path):
    displacements = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            atom_index = int(parts[0])
            displacement = [float(parts[1]), float(parts[2]), float(parts[3])]
            displacements[atom_index] = displacement
    return displacements

def apply_displacements(atomic_positions, displacements):
    modified_positions = []
    for i, line in enumerate(atomic_positions):
        parts = line.split()
        atom_index = i + 1
        if atom_index in displacements:
            displacement = displacements[atom_index]
            x = float(parts[1]) + displacement[0]
            y = float(parts[2]) + displacement[1]
            z = float(parts[3]) + displacement[2]
            modified_positions.append(f"{parts[0]} {x:.8f} {y:.8f} {z:.8f}\n")
        else:
            modified_positions.append(line)
    return modified_positions

def write_qe_input(file_path, lines, modified_positions, atomic_positions_start, atomic_positions_end):
    with open(file_path, 'w') as file:
        file.writelines(lines[:atomic_positions_start])
        file.writelines(modified_positions)
        file.writelines(lines[atomic_positions_end:])

def modify_qe_input(original_input_file, displacement_file, modified_input_file):
    lines, atomic_positions, start, end = read_qe_input(original_input_file)
    displacements = read_displacements(displacement_file)
    modified_positions = apply_displacements(atomic_positions, displacements)
    write_qe_input(modified_input_file, lines, modified_positions, start, end)

# Example usage:
for QE_input in QE_files_to_modify:
    original_file = f"../step_{current_step-1}/" + QE_input
    new_file = QE_input
    modify_qe_input(original_file, displacements_file, new_file)

print('Finished!')

