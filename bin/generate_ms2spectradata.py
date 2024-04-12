import json
import argparse

def generate_ms2data(input):
    ms2_data = {}
    with open(input, 'r') as input_file:
        add = 0
        peaks = []
        scan = 0
        charge = 0
        in_peaks = False
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                in_spectra = True
            elif in_spectra and line.startswith("END IONS"):
                if add:
                    ms2_data[scan] = {"Charge" : charge, "Peaks" : peaks}
                peaks = []
                add = False
                in_peaks = False
                in_spectra = False
            elif in_spectra and line.startswith("MSLEVEL"):
                add = (int(line.split("=")[1]) == 2)
            elif in_spectra and line.startswith("SCANS"):
                scan = int(line.split("=")[1])
            elif in_spectra and line.startswith("CHARGE") and add:
                charge = int(line.split("=")[1])
            elif in_spectra and line.startswith("Num peaks") and add:
                in_peaks = True
            elif in_spectra and in_peaks:
                data = line.split(" ")
                mass, intensity = float(data[0]), float(data[1])
                peaks.append(tuple([mass, intensity]))
            else:
                continue
    return ms2_data


def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('mgf')
    parser.add_argument('output_filename')
    args = parser.parse_args()

    spectraData = generate_ms2data(args.mgf)
    with open(args.output_filename, 'w') as file:
        json.dump(spectraData, file, indent = 2)
    
if __name__ == "__main__":
    main()

