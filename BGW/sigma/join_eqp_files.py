import os
import sys
import argparse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_eqp_file import read_eqp_dat_file

'''
Join multiple BerkeleyGW eqp.dat files (e.g. eqp1.dat from sigma.cplx.x)
that cover different, non-overlapping band ranges at the SAME k-points
into a single eqp_all.dat file.

Motivation: BerkeleyGW's sigma.inp takes one contiguous
band_index_min/band_index_max range. When the desired band range is too
expensive to run in a single job, split it into 2+ sigma runs over
sub-ranges (each producing its own eqp1.dat with the same k-point list),
then join them here into one file covering the full range.

Usage:
    python join_eqp_files.py -files eqp_part1.dat eqp_part2.dat eqp_part3.dat -out eqp_all.dat

Output format matches BerkeleyGW's own eqp.dat convention (one header
line "kx ky kz Nbnd_total" per k-point, followed by Nbnd_total lines of
"spin band_index E_DFT E_QP"), so eqp_all.dat can be used directly
wherever a normal eqp.dat/eqp1.dat is expected (e.g. sig2wan.x).
'''


def check_kpoints_match(kpoints_list, file_names, tol=1e-6):
    ref = kpoints_list[0]
    for i in range(1, len(kpoints_list)):
        kp = kpoints_list[i]
        if kp.shape != ref.shape:
            raise ValueError(
                f"'{file_names[i]}' has {kp.shape[0]} k-points, "
                f"expected {ref.shape[0]} (from '{file_names[0]}')"
            )
        if not np.allclose(kp, ref, atol=tol):
            bad = np.where(~np.all(np.isclose(kp, ref, atol=tol), axis=1))[0]
            raise ValueError(
                f"'{file_names[i]}'s k-points do not match '{file_names[0]}'s "
                f"(tol={tol}). First mismatching k-point index: {bad[0]}\n"
                f"  {file_names[0]}: {ref[bad[0]]}\n"
                f"  {file_names[i]}: {kp[bad[0]]}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Join multiple BerkeleyGW eqp.dat files covering "
                    "different band ranges (same k-points) into one eqp_all.dat"
    )
    parser.add_argument("-files", "--eqp_files", nargs="+", required=True,
                         help="eqp.dat files to join, any order (e.g. -files "
                              "eqp_part1.dat eqp_part2.dat eqp_part3.dat)")
    parser.add_argument("-out", "--out_file", default="eqp_all.dat",
                         help="Output file name (default: eqp_all.dat)")
    args = parser.parse_args()

    # 1. Load Eqp, Edft, band_indexes, Kpoints_list for each file
    all_data = []
    for f in args.eqp_files:
        print(f"Reading {f}")
        bands_dft, bands_qp, Kpoints, Nk, band_indexes = read_eqp_dat_file(f)
        band_indexes = band_indexes.astype(int)
        print(f"  {Nk} k-points, bands {band_indexes.min()}-{band_indexes.max()} "
              f"({len(band_indexes)} bands)")
        all_data.append(dict(file=f, bands_dft=bands_dft, bands_qp=bands_qp,
                              Kpoints=Kpoints, Nk=Nk, band_indexes=band_indexes))

    # 2. Certify k-points are the same across all files
    check_kpoints_match([d["Kpoints"] for d in all_data], args.eqp_files)
    Kpoints = all_data[0]["Kpoints"]
    Nk = all_data[0]["Nk"]
    print(f"\nAll {len(all_data)} files share the same {Nk} k-points -- OK")

    # 3. Build combined band index list in increasing order, dropping
    # repeated indices (keep the FIRST file's data, warn about it), and
    # warn about any gaps in the resulting range
    band_to_source = {}  # band_index -> (file_idx, local_row_idx)
    duplicates = {}
    for file_idx, d in enumerate(all_data):
        for local_idx, b in enumerate(d["band_indexes"]):
            if b in band_to_source:
                duplicates.setdefault(b, [band_to_source[b][0]]).append(file_idx)
            else:
                band_to_source[b] = (file_idx, local_idx)

    if duplicates:
        print("\nWARNING: duplicate band indices found across files "
              "(keeping the FIRST file listed for each, in -files order):")
        for b in sorted(duplicates):
            files_involved = duplicates[b]
            kept_idx = band_to_source[b][0]
            names = [all_data[i]["file"] for i in files_involved]
            print(f"  band {b}: present in {names} -- using data from "
                  f"'{all_data[kept_idx]['file']}'")

    band_indexes_total = np.array(sorted(band_to_source.keys()))

    full_range = np.arange(band_indexes_total.min(), band_indexes_total.max() + 1)
    missing = sorted(set(full_range.tolist()) - set(band_indexes_total.tolist()))
    if missing:
        print(f"\nWARNING: missing band indices in the combined range "
              f"{band_indexes_total.min()}-{band_indexes_total.max()}: {missing}")
    else:
        print(f"\nCombined band range {band_indexes_total.min()}-"
              f"{band_indexes_total.max()} is contiguous, no gaps -- OK")

    print(f"\nTotal unique bands to write: {len(band_indexes_total)}")

    # 4. Write eqp_all.dat, BerkeleyGW eqp.dat format
    Nbnd_total = len(band_indexes_total)
    with open(args.out_file, "w") as fout:
        for ik in range(Nk):
            kx, ky, kz = Kpoints[ik]
            fout.write(f" {kx:14.9f} {ky:14.9f} {kz:14.9f} {Nbnd_total:8d}\n")
            for b in band_indexes_total:
                file_idx, local_idx = band_to_source[int(b)]
                d = all_data[file_idx]
                e_dft = d["bands_dft"][local_idx, ik]
                e_qp = d["bands_qp"][local_idx, ik]
                fout.write(f"{1:8d}{int(b):8d}{e_dft:16.9f}{e_qp:16.9f}\n")

    print(f"\nWrote {args.out_file}")


if __name__ == "__main__":
    main()
