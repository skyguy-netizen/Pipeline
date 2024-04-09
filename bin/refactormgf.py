import argparse

def refactor_mgf(input, output):
    with open(input, 'r') as input_file, open(output, 'w') as output_file:
        spectra = []
        in_spectra = False
        feature_id = 0
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                in_spectra = True
                spectra.append(line)
            elif in_spectra and line.startswith("END IONS"):
                spectra.append(line)
                in_spectra = False
                output_file.write("\n".join(spectra))
                output_file.write("\n\n")
                feature_id += 1
                spectra = []
            elif in_spectra and line.startswith("FEATURE_ID") or line.startswith("SCANS"):
                updated_fid = line.split("=")
                updated_fid[-1] = str(feature_id)
                spectra.append("=".join(updated_fid))
            elif in_spectra:
                spectra.append(line)
    return

def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')

    args = parser.parse_args()
    refactor_mgf(args.input_filename, args.output_filename)

if __name__ == "__main__":
    main()
