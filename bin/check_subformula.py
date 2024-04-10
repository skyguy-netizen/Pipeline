import argparse
import json
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as c
import re
from collections import defaultdict

def get_elements(smile):
    try:
        formula = get_formula(get_molecule(smile))
    except:
        return defaultdict(int)
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    elements = defaultdict(int)
    for element, count in elements:
        elements[element] = count
    return elements

def get_formula(mol):
    return c.CalcMolFormula(mol)

def get_molecule(smile):
    return Chem.MolFromSmiles(smile)

def isSubFormula(predElems, originalElems):
    for k in predElems.keys():
        if predElems[k] > originalElems[k]:
            return False
    return True

def checkSubFormula(tree : dict):
    for smil_add, scans in tree.items():
        smile = json.loads(smil_add)[0]
        for ms2, ms3_scans in scans.items():
            newms3 = []
            for scan in ms3_scans:
                if not scan['Substructure?']:
                    pred_elements = get_elements(scan['PredSmile'])
                    orig_elements = get_elements(smile)
                    # print(type)
                    scan['Subformula'] = isSubFormula(pred_elements, orig_elements)
                newms3.append(scan)
            tree[smil_add][ms2] = newms3
    return


def main():
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('tree_file')
    parser.add_argument('output_filename')
    args = parser.parse_args()

    tree = json.load(open(args.tree_file, 'r'))
    checkSubFormula(tree)
    with open(args.output_filename, 'w') as file:
        json.dump(tree, file, indent = 2)

if __name__ == "__main__":
    main()