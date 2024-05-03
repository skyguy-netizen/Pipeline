import argparse
from rdkit import Chem
import json

def is_same(smile1, smile2):
    if smile2 == float('nan') :
        return True
    
    # Convert SMILES strings to RDKit molecules
    try:
        molecule1 = Chem.MolFromSmiles(smile1)
        molecule2 = Chem.MolFromSmiles(smile2)
    except:
        return True

    if not molecule1 and not molecule2:
        return True
    elif (not molecule2):
        return True
    elif (not molecule1):
        return False
    else:
    # try:
    # # check if they are the same
        return (molecule1.HasSubstructMatch(molecule2) and molecule2.HasSubstructMatch(molecule1))
    
def removeDuplicates(tree):
    for mol, scans in tree.items():
        for ms2, ms3_scans in scans.items():
            new_scans = []
            for scan in ms3_scans:
                add = True
                for j in new_scans:
                    # print()
                    # print()
                    # print(type(scan['PredSmile']))
                    # print(j)
                    if is_same(j['PredSmile'], scan['PredSmile']):
                        j['Substructure?'] = j['Substructure?'] or scan['Substructure?']
                        add = False
                        break
                if add:
                    if scan['PredSmile'] is not None:
                        new_scans.append(scan)
                else:
                    continue
            tree[mol][ms2] = new_scans

    

def main():    
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('tree_file')
    parser.add_argument('output_filename')
    args = parser.parse_args()

    tree = json.load(open(args.tree_file, 'r'))
    removeDuplicates(tree)

    with open(args.output_filename, 'w') as file:
        json.dump(tree, file, indent = 2)  

if __name__ == "__main__":
    main()