import json
import argparse

def merge(tree, spectradata):
    print("In merge")
    for smil_add, scans in tree.items():
        for ms2, ms3_scans in scans.items():
            if ms2 == "Not found scan":
                continue
            spectradata[ms2]['MS3 Scans'] = ms3_scans
            tree[smil_add][ms2] = spectradata[ms2]
    return

def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('final_tree')
    parser.add_argument('ms2_spectradata')
    args = parser.parse_args()

    tree = json.load(open(args.final_tree, 'r'))
    spectradata = json.load(open(args.ms2_spectradata, 'r'))
    merge(tree, spectradata)

    with open(args.final_tree, 'w') as file:
        json.dump(tree, file, indent = 1)

if __name__ == "__main__":
    main()