import argparse
import pandas as pd
import json
import requests
from rdkit import Chem

url = "https://gnps2.org/result?task=121341ce6b534fff8bb7898720923a7c&viewname=summarylist&json"

def getscan(smile, adduct, energy, ms2scans):
    filtered_row = ms2scans.loc[(ms2scans["Smile"] == smile) & (ms2scans["Adduct"] == adduct) & (ms2scans["Energies"] == energy)]
    if not filtered_row.empty:
        scans_list = filtered_row.iloc[0]["Scan"]
        return scans_list
    else:
        return "No scan found"

def is_subgraph(original_smiles, predicted_smiles):
    # Convert SMILES strings to RDKit molecules
    original_mol = Chem.MolFromSmiles(original_smiles)
    predicted_mol = Chem.MolFromSmiles(predicted_smiles)

    # Check if the predicted molecule is a subgraph of the original molecule
    return original_mol.HasSubstructMatch(predicted_mol)

def generate_ms3data(ms3input, full_dict):
    with open(ms3input, 'r') as input_file:
        spectra = {}
        in_spectra = False
        ms2_peak = ""
        smile = ""
        adduct = ""
        scan = 0
        col_energies = []
        masses = []
        for line in input_file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                in_spectra = True
            elif in_spectra and line.startswith("END IONS"):
                in_spectra = False
                try:
                    spectra[smile]["Actual"][scan] = ms2_peak
                except:
                    spectra[smile] = {"Actual" : {scan : ms2_peak}}
                full_dict[scan] = [ms2_peak, smile, col_energies, adduct, masses]
            elif in_spectra and line.startswith("SMILES"):
                smile = line[7:]
            elif in_spectra and line.startswith("PEPMASS"):
                ms2_peak = line.split("=")[1]
            elif in_spectra and line.startswith("MSn_collision_energies"):
                ms1_energy = float(line.split("=")[1].split(",")[0].strip()[1:])
                ms2_energy = float(line.split("=")[1].split(",")[1].strip()[0:-1])
                col_energies = [ms1_energy, ms2_energy]
            elif in_spectra and line.startswith("MSn_precursor_mzs"):
                ms1_mass = float(line.split("=")[1].split(",")[0].strip()[1:])
                ms2_mass = float(line.split("=")[1].split(",")[1].strip()[0:-1])
                masses = [ms1_mass, ms2_mass]
            elif in_spectra and line.startswith("ADDUCT"):
                adduct = line[7:]
            elif in_spectra and line.startswith("SCANS"):
                scan = line.split("=")[1]
    return spectra


def combine_results(workflowdata, full_dict, refactored_data):
    for data_list in workflowdata:
        scan_num = str(data_list['#Scan#'])
        pred_smile = data_list["Smiles"]
        
        smile = full_dict[scan_num][1]
        energies = full_dict[scan_num][2]
        adduct = full_dict[scan_num][3]
        masses = full_dict[scan_num][4]

        #dump tuple to json string
        map_id = json.dumps((smile, adduct))

        try:
            refactored_data[map_id].append({"Scan" : scan_num, "PredSmile" : pred_smile, "Energies" : energies, "PrecMasses" : masses})
        except:
                refactored_data[map_id] = [{"Scan" : scan_num, "PredSmile" : pred_smile, "Energies" : energies, "PrecMasses" : masses}]

    return
    
def check_substructure(refactored_data, ms2scans, output):
    final_results_grouped = {}
    for og_map_id, ms3_spectras in refactored_data.items():
        smile_adduct = json.loads(og_map_id)
        for spectra in ms3_spectras:
            check = False
            if(spectra['PredSmile'] is not None):
                try:
                    check = is_subgraph(smile_adduct[0], spectra['PredSmile'])
                except:
                    pass
            spectra['Substructure?']  = check


            try:
                ms2_scan = int(getscan(smile_adduct[0], smile_adduct[1], spectra["Energies"][0], ms2scans)[0])
            except:

                ms2_scan = "Not found scan"
            
            try:
                final_results_grouped[og_map_id][ms2_scan].append(spectra)
            except:
                try:
                    final_results_grouped[og_map_id][ms2_scan] = [spectra]
                except:
                    final_results_grouped[og_map_id] = {ms2_scan : [spectra]}
    
    with open(output, 'w') as f:
        json.dump(final_results_grouped, f, indent=2)


def main():
    
    full_dict = {}  
    refactored_data = {}
    
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('ms3input')
    parser.add_argument('ms2scans')
    parser.add_argument('workflow_results')
    parser.add_argument('output_filename')

    args = parser.parse_args()

    ms2scans = pd.read_json(args.ms2scans, orient = 'records')
    workflow_data = pd.read_table(args.workflow_results, sep = '\t').to_dict(orient = 'records')
    generate_ms3data(args.ms3input, full_dict)
    combine_results(workflow_data, full_dict, refactored_data)
    check_substructure(refactored_data, ms2scans, args.output_filename)



if __name__ == "__main__":
    main()
