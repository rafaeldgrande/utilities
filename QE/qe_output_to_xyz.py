import argparse
from ase.io import read, write

def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert Quantum ESPRESSO output to XYZ format using ASE."
    )
    parser.add_argument("input", help="QE output file (e.g. relax.out)")
    parser.add_argument("output", help="XYZ output file (e.g. output.xyz)")
    parser.add_argument(
        "--index", default="-1",
        help="Frame index or slice (default: -1 = last frame). Use ':' for all frames."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    index = args.index if args.index == ":" else int(args.index)
    atoms = read(args.input, format="espresso-out", index=index)
    write(args.output, atoms)
    print(f"Written: {args.input} → {args.output}  (index={args.index})")

if __name__ == "__main__":
    main()
