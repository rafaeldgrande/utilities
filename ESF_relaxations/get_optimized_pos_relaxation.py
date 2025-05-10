
import sys

def extract_last_atomic_positions(file_path):
    """
    Extracts the last occurrence of atomic positions from a file containing 
    multiple "ATOMIC_POSITIONS" sections.

    The function reads the file line by line, identifies all occurrences of 
    the "ATOMIC_POSITIONS" section, and returns the atomic positions from 
    the last occurrence.

    Args:
        file_path (str): The path to the file containing the atomic positions.

    Returns:
        list of str: A list of strings representing the atomic positions from 
        the last "ATOMIC_POSITIONS" section. Each string corresponds to a line 
        in the section.
    """
    # Read the file and search for the last ATOMIC_POSITIONS section
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atomic_positions = []
    temp_positions = []
    recording = False

    # Loop through the lines to find all occurrences of the ATOMIC_POSITIONS section
    for line in lines:
        if line.strip().startswith('ATOMIC_POSITIONS'):
            recording = True
            temp_positions = []  # Reset for new atomic positions
            continue

        if recording:
            if line.strip() == '':  # Stop recording when an empty line is reached
                recording = False
                atomic_positions = temp_positions  # Keep the last set
            else:
                temp_positions.append(line.strip())

    return atomic_positions

def write_atomic_positions(atomic_positions, output_file):
    """
    Writes a list of atomic positions to a specified output file in the format
    required for atomic position data.

    The output file will include a header "ATOMIC_POSITIONS (angstrom)" followed
    by each atomic position on a new line.

    Args:
        atomic_positions (list of str): A list of strings where each string 
            represents an atomic position.
        output_file (str): The path to the output file where the atomic positions 
            will be written.

    Returns:
        None
    """
    # Write the atomic positions to a new file
    with open(output_file, 'w') as file:
        file.write("ATOMIC_POSITIONS (angstrom)\n")
        for position in atomic_positions:
            file.write(position + "\n")

if __name__ == "__main__":
    # Check if the input file path is provided
    if len(sys.argv) < 3:
        print("Usage: python extract_positions.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Extract the last atomic positions
    atomic_positions = extract_last_atomic_positions(input_file)

    # Write the atomic positions to the output file
    write_atomic_positions(atomic_positions, output_file)

    print(f"Optimized atomic positions have been written to {output_file}")

