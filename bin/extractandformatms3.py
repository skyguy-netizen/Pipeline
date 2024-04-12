import argparse

def extractms3(input, intermediate):
    with open(input, 'r') as input_file, open(intermediate, 'w') as output_file:
        spectra = []
        add = False
        in_spectra = False
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                add = False
                in_spectra = True
                spectra.append(line)
            elif in_spectra and line.startswith("END IONS"):
                spectra.append(line)
                in_spectra = False
                if add:
                    output_file.write("\n".join(spectra))
                    output_file.write("\n\n")
                spectra = []
                add = False
            elif in_spectra and line.startswith("MSLEVEL"):
                ms_level = line.split("=")[1]
                if (ms_level == "3"):
                    add = True
                spectra.append(line)
            elif in_spectra:
                spectra.append(line)
    return

def formatms3(intermediate, output):
    with open(intermediate, 'r') as input_file, open(output, 'w') as output_file:
        spectra = []
        add = False
        ms2_spectra = 0.0
        in_spectra = False
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                add = False
                in_spectra = True
                spectra.append(line)
            elif in_spectra and line.startswith("END IONS"):
                spectra.append(line)
                in_spectra = False
                if add:
                    output_file.write("\n".join(spectra))
                    output_file.write("\n\n")
                spectra = []
                add = False
            elif in_spectra and line.startswith("PEPMASS"):
                continue
            elif in_spectra and line.startswith("MSn_precursor_mzs"):
                ms2_spectra = float(line.split("=")[1].split(",")[1].strip()[0:-1])
                spectra.append(line)
                spectra.insert(12, "PEPMASS=" + str(ms2_spectra))
            elif in_spectra and line.startswith("MSLEVEL"):
                ms_level = line.split("=")[1]
                if (ms_level == "3"):
                    add = True
                spectra.append(line)
            elif in_spectra:
                spectra.append(line)
    return None

def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')

    args = parser.parse_args()
    extractms3(args.input_filename, 'intermediate.json')
    formatms3('intermediate.json', args.output_filename)

if __name__ == "__main__":
    main()
