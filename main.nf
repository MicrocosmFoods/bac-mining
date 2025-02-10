#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=16
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL GENOMES FOR PEPTIDES AND BGCS, AND PERFORM
FUNCTIONAL ANNOTATION.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH AND KOFAM DATABASES AND PROVIDE THE PATHS
WITH --antismash_db AND --kofam_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING 
DATABASES AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
antismash_db                    : $params.antismash_db
kofam_db                        : $params.kofam_db
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels
// genome_fastas tuple with genome name and fasta filepath
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

antismash_db_ch = channel.fromPath(params.antismash_db)
kofam_db_ch = channel.fromPath(params.kofam_db)

// workflow steps
workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[0] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    genome_stb_tsv = make_genome_stb.out.stb_tsv
    
    // get small ORF predictions with smorfinder and combine into a single FASTA
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file.collect()
    smorfinder_tsvs = smorfinder.out.tsv_file.collect()
    combine_smorf_proteins(smorf_proteins)
    all_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // predict ORFs with pyrodigal
    pyrodigal(genome_fastas)
    predicted_orfs_gbks = pyrodigal.out.predicted_orfs_gbk
    predicted_orfs_proteins = pyrodigal.out.predicted_orfs_faa

    // predict cleavage peptides with deeppeptide, extract sequences from json, and combine into a single FASTA
    predict_cleavage_peptides(predicted_orfs_proteins)
    cleavage_peptides_json = predict_cleavage_peptides.out.cleavage_peptides_json
    cleavage_input_ch = cleavage_peptides_json.join(predicted_orfs_proteins, by: 0)
    extract_cleavage_peptides_json(cleavage_input_ch)
    all_cleavage_peptides_fastas = extract_cleavage_peptides_json.out.cleavage_peptides_fasta.collect()
    deeppeptide_tsvs = extract_cleavage_peptides_json.out.cleavage_peptides_tsv.collect()
    combine_cleavage_peptides(all_cleavage_peptides_fastas)
    all_cleavage_peptides = combine_cleavage_peptides.out.combined_cleavage_peptides

    // antismash predictions
    antismash_input_ch = predicted_orfs_gbks.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results
    all_antismash_gbk_files = antismash_gbk_files.map{ it[1] }.collect()

    // extract antismash gbks into summary tsv, ripp peptides
    extract_gbks(all_antismash_gbk_files, genome_stb_tsv)
    antismash_summary_tsv = extract_gbks.out.bgc_summary_tsv
    antismash_peptides_tsv = extract_gbks.out.bgc_peptides_tsv

    // summarize peptide and BGC counts per genome
    summarize_molecule_counts(smorfinder_tsvs, deeppeptide_tsvs, antismash_summary_tsv, antismash_peptides_tsv)

    // run kofamscan annotations on all predicted proteins
    kofamscan_annotation_ch = predicted_orfs_proteins.combine(kofam_db_ch)
    kofamscan_annotation(kofamscan_annotation_ch)
    all_kofamscan_tsvs = kofamscan_annotation.out.kofamscan_tsv.collect()

    // combine kofamscan results
    combine_kofamscan_results(all_kofamscan_tsvs)
}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv

    script:
    """
    python ${baseDir}/bin/generate_genome_stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    errorStrategy 'ignore'
    // rarely some genomes will fail for no discernible reason, skip over these

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:v0.2"

    input:
    tuple path(fasta), val(genome_name)

    output:
    path("*.gff"), emit: gff_file
    path("*.faa"), emit: faa_file
    path("*.ffn"), emit: ffn_file
    path("*.tsv"), emit: tsv_file

    script:
    """
    smorf single ${fasta} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff
    ln -s ${genome_name}/${genome_name}.faa
    ln -s ${genome_name}/${genome_name}.ffn
    ln -s ${genome_name}/${genome_name}_smorfinder.tsv
    """

}

process combine_smorf_proteins {
    tag "combine_smorf_proteins"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = '10 GB'
    cpus = 1
    
    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(smorf_proteins)

    output: 
    path("*.fasta"), emit: combined_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} combined_smorf_proteins.fasta
    """
    
}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/pyrodigal:3.4.1--py310h4b81fae_0"

    input:
    tuple path(fasta), val(genome_name)

    output:
    tuple val(genome_name), path("*.fna"), emit: predicted_orfs_fna
    tuple val(genome_name), path("*.faa"), emit: predicted_orfs_faa
    tuple val(genome_name), path("*.gbk"), emit: predicted_orfs_gbk

    script:
    """
    pyrodigal \\
        -i ${fasta} \\
        -f "gbk" \\
        -o "${genome_name}.gbk" \\
        -d ${genome_name}.fna \\
        -a ${genome_name}.faa \\
        --no-stop-codon
    """
}

process predict_cleavage_peptides {
    tag "${genome_name}_predict_cleavage_peptides"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    accelerator 1, type: 'nvidia-t4'
    cpus = 8

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/deeppeptide:v0.5"

    input:
    tuple val(genome_name), path(predicted_orfs_proteins)

    output:
    tuple val(genome_name), path("*.json"), emit: cleavage_peptides_json

    script:
    """
    WORKDIR=\$PWD
    cp ${predicted_orfs_proteins} /app/DeepPeptide/predictor
    cd /app/DeepPeptide/predictor

    python3 predict.py --fastafile ${predicted_orfs_proteins.getName()} --output_dir ${genome_name} --output_fmt json --batch_size 500
    
    cp ${genome_name}/*.json \$WORKDIR/
    """
}

process extract_cleavage_peptides_json {
    tag "${genome_name}_extract_cleavage_peptides_json"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    tuple val(genome_name), path(deeppeptide_json), path(protein_faa)

    output:
    path("*_parent_proteins.faa"), emit: parent_proteins_faa
    path("*_peptides.faa"), emit: cleavage_peptides_fasta
    path("*.tsv"), emit: cleavage_peptides_tsv

    script:
    """
        python ${baseDir}/bin/extract_cleavage_peptides_json.py \
            --json_file ${deeppeptide_json} \
            --protein_fasta_file ${protein_faa} \
            --proteins_output_file ${genome_name}_parent_proteins.faa \
            --protein_peptides_output_file ${genome_name}_peptides.faa \
            --predictions_output_file ${genome_name}_deeppeptide.tsv
    """
}

process combine_cleavage_peptides {
    tag "combine_cleavage_peptides"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(cleavage_peptides)

    output:
    path("*.fasta"), emit: combined_cleavage_peptides

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${cleavage_peptides.join(' ')} combined_cleavage_peptides.fasta
    """
}


process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    memory = "20 GB"
    cpus = 6

    container "public.ecr.aws/biocontainers/antismash-lite:7.1.0--pyhdfd78af_0"

    input:
    tuple val(genome_name), path(gbk_file), path(databases)

    output: 
    tuple val(genome_name), path("${genome_name}/*.json") , emit: json_results
    tuple val(genome_name), path("${genome_name}/*.log") , emit: log
    tuple val(genome_name), path("${genome_name}/*region*.gbk") , optional: true, emit: gbk_results

    script: 
    """
    antismash \\
        -c $task.cpus \\
        --output-dir ${genome_name} \\
        --output-basename ${genome_name} \\
        --logfile ${genome_name}/${genome_name}.log \\
        --databases $databases \\
        --genefinding-tool none \\
        ${gbk_file}
    """
}

process extract_gbks {
    tag "extract_gbks"
    publishDir "${params.outdir}/main_results/bgc_info", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input: 
    path(gbk_files)
    path(genome_stb_tsv)

    output:
    path("antismash_summary.tsv"), emit: bgc_summary_tsv
    path("antismash_peptides.tsv"), emit: bgc_peptides_tsv, optional: true
    path("antismash_peptides.fasta"), emit: bgc_peptides_fasta, optional: true

    script:
    """
    python ${baseDir}/bin/extract_antismash_gbks.py ${gbk_files.join(' ')} ${genome_stb_tsv} antismash_summary.tsv antismash_peptides.tsv antismash_peptides.fasta
    """

}

process summarize_molecule_counts {
    tag "summarize_molecule_counts"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(smorfinder_tsvs)
    path(deeppeptide_tsvs)
    path(antismash_summary_tsv)
    path(antismash_peptides_tsv)

    output:
    path("all_molecule_counts.tsv"), emit: all_molecule_counts_tsv
    path("all_smorfinder_results.tsv"), emit: all_smorfinder_results_tsv
    path("all_deeppeptide_results.tsv"), emit: all_deeppeptide_results_tsv

    script:
    """
    python3 ${baseDir}/bin/process_molecule_counts.py \\
    --smorfinder-tsvs ${smorfinder_tsvs.join(' ')} \\
    --deeppeptide-tsvs ${deeppeptide_tsvs.join(' ')} \\
    --bgc-summary ${antismash_summary_tsv} \\
    --antismash-peptides ${antismash_peptides_tsv} \\
    --output-counts all_molecule_counts.tsv \\
    --output-smorfinder all_smorfinder_results.tsv \\
    --output-deeppeptide all_deeppeptide_results.tsv
    """
}

process kofamscan_annotation {
    tag "${genome_name}_kofam_scan_annotation"
    publishDir "${params.outdir}/kofam_scan_annotation", mode: 'copy'

    memory = "25 GB"
    cpus = 12

    container "public.ecr.aws/biocontainers/kofamscan:1.0.0--0"

    input:
    tuple val(genome_name), path(faa_file), path(kegg_db_dir)

    output:
    path("*.tsv"), emit: kofamscan_tsv

    script:
    """
    exec_annotation --format detail --ko-list ${kegg_db_dir}/ko_list --profile ${kegg_db_dir}/profiles --cpu ${task.cpus} -o ${genome_name}_kofamscan_annotations.tsv ${faa_file}
    """
}

process combine_kofamscan_results {
    tag "combine_kofamscan_results"
    publishDir "${params.outdir}/main_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/polars:0.12.5"

    input:
    path(kofamscan_tsvs)

    output:
    path("combined_kofamscan_results.tsv"), emit: combined_kofamscan_results_tsv

    script:
    """
    python3 ${baseDir}/bin/combine_kofamscan_results.py \\
    --input_files ${kofamscan_tsvs.join(' ')} \\
    --output combined_kofamscan_results.tsv
    """
}