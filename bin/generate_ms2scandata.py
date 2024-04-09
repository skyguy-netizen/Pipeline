import argparse
import pandas as pd

def generate_ms2data(input):
    ms2_data = []
    with open(input, 'r') as input_file:
        add = 0
        spectra = {}
        spectra["SpecType"] = None
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                in_spectra = True
            elif in_spectra and line.startswith("END IONS"):
                if add:
                    ms2_data.append(spectra)
                spectra = {}
                spectra["SpecType"] = None
                add = False
            elif in_spectra and line.startswith("SMILES"):
                spectra["Smile"] = line[7:]
            elif in_spectra and line.startswith("MSLEVEL"):
                ms = int(line.split("=")[1])
                add = (ms == 2)
            elif in_spectra and line.startswith("ADDUCT"):
                spectra["Adduct"] = line[7:]
            elif in_spectra and line.startswith("FEATURE_ID"):
                spectra["Scan"] = int(line.split("=")[-1])
            elif in_spectra and line.startswith("Collision energy") and add:
                spectra["Energies"] = float(line.split("=")[-1])
            elif in_spectra and line.startswith("SPECTYPE") and add:
                spectra["SpecType"] = line[9:]
            elif in_spectra and line.startswith("CHARGE") and add:
                spectra["precursor_charge"] = float(line.split("=")[-1])
            elif in_spectra and line.startswith("PEPMASS") and add:
                spectra["precursor_mass"] = float(line.split("=")[-1])
            else:
                continue
    return ms2_data


def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')
    args = parser.parse_args()

    df = pd.DataFrame(generate_ms2data(args.input_filename))
    df = df.sort_values(by = "Energies").reset_index(drop = True)
    ms2scans = df.groupby(["Smile", "Adduct", "Energies"])["Scan"].apply(list).reset_index()
    ms2scans.to_csv(args.output_filename)


if __name__ == "__main__":
    main()