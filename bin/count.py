import argparse
import json
from collections import defaultdict



def count(input):
    counts_ms3 = defaultdict(int)
    counts_ms2 = defaultdict(int)
    smiles = set()
    with open(input, 'r') as input_file:
        smile = ""
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                in_spectra = True
            elif in_spectra and line.startswith("END IONS"):
                in_spectra = False
            elif in_spectra and line.startswith("SMILES"):
                smile = line[7:]
                smiles.add(smile)
            elif in_spectra and line.startswith("MSLEVEL"):
                ms = int(line.split("=")[1])
                if ms == 3:
                    counts_ms3[smile] += 1
                elif ms == 2:
                    counts_ms2[smile] += 1
                else:
                    continue
            else:
                continue
    return (counts_ms3, counts_ms2, smiles)



def combine_dicts(ms3, ms2, output, smiles):
    combined = {}
    for smile in smiles:
        combined[smile] = {"Total MS3 cnts" : ms3[smile], "Total MS2 cnts" : ms2[smile]}

    with open(output, 'w') as f:
        json.dump(combined, f, indent = 2)

    return



def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')

    args = parser.parse_args()
    (ms3, ms2, smiles) = count(args.input_filename)
    combine_dicts(ms3, ms2, args.output_filename, smiles)


if __name__ == "__main__":
    main()