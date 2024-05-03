#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { LibrarySearch } from "$baseDir/LibrarySearch_Workflow/nf_workflow.nf"

params.input_merge = "$baseDir/data/*.mgf"
params.inputlibraries = "$baseDir/data/libraries"
params.inputspectra = "$baseDir/data/spectra"

TOOL_FOLDER = "$baseDir/bin"
INPUT_FOLDER = "$baseDir/data"
OUTPUT_FOLDER = "$baseDir/nf_output"

process processOriginalMGF {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input

    output:
    file '*.mgf'

    """
    python $TOOL_FOLDER/refactormgf.py $input ${input.baseName + '-refactored'}.mgf
    """
}


process extractMS3data {
    publishDir "./data/spectra", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input

    output:
    path "ms3.mgf", emit: ms3file

    """
    python $TOOL_FOLDER/extractandformatms3.py $input ms3.mgf 
    """
}

process counts {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input

    output:
    file "counts.json"

    """
    python $TOOL_FOLDER/count.py $input counts.json
    """
}

process generate_ms2scandata {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file original

    output:
    file 'ms2scandata.json'

    """
    python $TOOL_FOLDER/generate_ms2scandata.py $original ms2scandata.json
    """
}

process combineWorkflowResults {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"
    
    input:
    file ms2scandata
    file ms3
    file workflowResults

    output:
    file 'combined_data.json'

    """
    python $TOOL_FOLDER/combineResults.py $ms2scandata $ms3  $workflowResults combined_data.json
    """ 
}

process removeDuplicates {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"
    
    input:
    file tree

    output:
    file 'unique_tree.json'

    """
    python $TOOL_FOLDER/removeDuplicates.py $tree unique_tree.json
    """
}

process checkSubFormula {
    publishDir "./nf_output", mode: 'copy'  

    conda "$TOOL_FOLDER/conda_env.yml"
    
    input:
    file tree

    output:
    file 'final_tree.json'

    """
    python $TOOL_FOLDER/check_subformula.py $tree final_tree.json
    """
}

process generateSpectraData {
    // publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file mgf

    output:
    file 'ms2spectradata.json'

    """
    python $TOOL_FOLDER/generate_ms2spectradata.py $mgf ms2spectradata.json
    """
}

process addSpectraData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file tree
    file spectradata

    output:
    file tree

    """
    python $TOOL_FOLDER/addms2spectradata.py $tree $spectradata
    """
}

// process runLibrarySearch {
//     publishDir "./nf_output", mode: 'copy'

//     conda "$TOOL_FOLDER/conda_env.yml"

//     input:
//     val ready

//     """
//     cd $baseDir/LibrarySearch_Workflow && nextflow run nf_workflow.nf -c nextflow.config -Dcapsule.log=verbose
//     """
// }

process dummy {
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val path

    """
    echo $path
    """
}

workflow {
    data = Channel.fromPath(params.input_merge, checkIfExists:true)
    results = processOriginalMGF(data)
    ms3 = extractMS3data(results.collate(1))    
    LibrarySearch(ms3.collate(1))
    //collect and flatten
    // dummy(extractMS3data.out.spectraDir)
    // LibrarySearch(extractMS3data.out.spectraFolder)
    count = counts(results.collate(1))
    workflow = generate_ms2scandata(results.collate(1))
    spectradata = generateSpectraData(results.collate(1))
    combined_data = combineWorkflowResults(ms3.collate(1), workflow.collate(1), workflowResults.collate(1))
    uniq_data = removeDuplicates(combined_data.collate(1))
    tree = checkSubFormula(uniq_data.collate(1))
    final_tree = addSpectraData(tree.collate(1), spectradata.collate(1))
}
