
import sys

def read_atomic_positions_from_qe_input(file_path):
    """Reads the atomic positions section from a Quantum Espresso input file."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atomic_positions = []
    recording = False

    # Loop through the lines to find the ATOMIC_POSITIONS section
    for line in lines:
        if line.strip().startswith('ATOMIC_POSITIONS'):
            recording = True
            continue

        if recording:
            if line.strip() == '':  # Stop recording at an empty line
                break
            atomic_positions.append(line.strip())

    return atomic_positions

def replace_atomic_positions_in_qe_input(input_file, new_positions, output_file):
    """Replaces the atomic positions in a Quantum Espresso input file."""
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Find where the atomic positions start and stop
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('ATOMIC_POSITIONS'):
            start_idx = i
            break

    if start_idx is not None:
        for i in range(start_idx + 1, len(lines)):
            if lines[i].strip() == '':  # Stop when encountering an empty line
                end_idx = i
                break

    # If we found the atomic positions section, replace it
    if start_idx is not None and end_idx is not None:
        # Keep everything before and after the atomic positions section
        new_lines = lines[:start_idx + 1] + [f"{pos}\n" for pos in new_positions] + lines[end_idx:]

        # Write the new content to the output file
        with open(output_file, 'w') as file:
            file.writelines(new_lines)
        print(f"Atomic positions successfully updated and written to {output_file}")
    else:
        print("Error: ATOMIC_POSITIONS section not found in the input file.")

def main(qe_input_file, atomic_positions_file, output_file):
    # Read new atomic positions from the second file
    new_atomic_positions = read_atomic_positions_from_qe_input(atomic_positions_file)

    # Replace the atomic positions in the QE input file and write to a new file
    replace_atomic_positions_in_qe_input(qe_input_file, new_atomic_positions, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python replace_atomic_positions.py <qe_input_file> <atomic_positions_file> <output_file>")
        sys.exit(1)

    qe_input_file = sys.argv[1]
    atomic_positions_file = sys.argv[2]
    output_file = sys.argv[3]

    main(qe_input_file, atomic_positions_file, output_file)

