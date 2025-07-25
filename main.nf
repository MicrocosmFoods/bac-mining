#! /usr/bin/env nextflow

// Description
// Mine bacterial genomes for different peptide types, BGCs, and perform functional annotation
// Note that for steps that process individual genomes, such as smorfinder, deeppeptide, antismash, etc. the error strategy is set to 'ignore'
// This is because individual genomes will sometimes fail for formatting reasons or not discernible reason, and we don't want these to halt the entire workflow

nextflow.enable.dsl=2

def date = new java.util.Date().format('yyyy-MM-dd')
params.outdir = "${date}-bacmagmining-results"
params.threads=16
params.functional_annotation = false
params.smorfinder_mode = 'single'

log.info """\

MINE BACTERIAL GENOMES FOR PEPTIDES AND BGCS, AND PERFORM
FUNCTIONAL ANNOTATION.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH AND KOFAM DATABASES AND PROVIDE THE PATHS
WITH --antismash_db AND --kofam_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING 
DATABASES AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
genome_list                     : $params.genome_list
antismash_db                    : $params.antismash_db
kofam_db                        : $params.kofam_db
outdir                          : $params.outdir
threads                         : $params.threads
functional_annotation           : $params.functional_annotation
smorfinder_mode                 : $params.smorfinder_mode
"""

// define channels
// genome_fastas tuple with genome name and fasta filepath
// if genome_list is provided, filter the genomes based on the list
if (params.genome_list) {
    genome_names = file(params.genome_list).readLines().collect { it.trim() }.toSet()
    genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
        .map { file -> 
            def baseName = file.getBaseName()
            return [baseName, file]
        }
        .filter { file, baseName -> genome_names.contains(baseName)}
} else {
    genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
        .map { file -> 
            def baseName = file.getBaseName()
            return [baseName, file]
        }
}

// Validate smorfinder_mode parameter
if (!['pre_called', 'single'].contains(params.smorfinder_mode)) {
    log.error "Invalid smorfinder_mode: ${params.smorfinder_mode}. Must be either 'pre_called' or 'single'"
    System.exit(1)
}

antismash_db_ch = channel.fromPath(params.antismash_db)
kofam_db_ch = channel.fromPath(params.kofam_db)

// workflow steps
workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[1] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    genome_stb_tsv = make_genome_stb.out.stb_tsv

    // predict ORFs with pyrodigal and save output files
    pyrodigal(genome_fastas)
    predicted_orfs_gbks = pyrodigal.out.predicted_orfs_gbk
    predicted_orfs_faa = pyrodigal.out.predicted_orfs_faa
    predicted_orfs_ffn = pyrodigal.out.predicted_orfs_ffn
    
    // convert .gbk to .gff for smorfinder input
    convert_gbk_to_gff(predicted_orfs_gbks)
    predicted_orfs_gff = convert_gbk_to_gff.out.predicted_orfs_gff
    
    // run smorfinder based on user preference
    if (params.smorfinder_mode == 'pre_called') {
        // call smorfinder with pre-called genes from main pyrodigal run
        smorfinder_input = genome_fastas.join(predicted_orfs_gff, by: 0).join(predicted_orfs_ffn, by: 0).join(predicted_orfs_faa, by: 0)
        smorfinder_pre_called(smorfinder_input)
        smorf_proteins = smorfinder_pre_called.out.smorf_faa.collect()
        smorfinder_tsvs = smorfinder_pre_called.out.smorf_tsv.collect()
    } else if (params.smorfinder_mode == 'single') {
        // call smorfinder directly on genomes (single mode)
        smorfinder_single(genome_fastas)
        smorf_proteins = smorfinder_single.out.smorf_faa.collect()
        smorfinder_tsvs = smorfinder_single.out.smorf_tsv.collect()
    }

    // combine smorfinder results and proteins
    combine_smorf_proteins(smorf_proteins)
    all_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // filter to proteins less than 50 AAs for inputting to deeppeptide
    filter_small_proteins(predicted_orfs_faa)
    filtered_proteins = filter_small_proteins.out.filtered_proteins

    // predict cleavage peptides with deeppeptide, extract sequences from json, and combine into a single FASTA
    predict_cleavage_peptides(filtered_proteins)
    cleavage_peptides_json = predict_cleavage_peptides.out.cleavage_peptides_json
    cleavage_input_ch = cleavage_peptides_json.join(predicted_orfs_faa, by: 0)
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
    summarize_molecule_counts(smorfinder_tsvs, deeppeptide_tsvs, antismash_summary_tsv)

    // Only run functional annotation if enabled
    if (params.functional_annotation) {
        // run kofamscan annotations on all predicted proteins
        kofamscan_annotation_ch = predicted_orfs_faa.combine(kofam_db_ch)
        kofamscan_annotation(kofamscan_annotation_ch)
        all_kofamscan_tsvs = kofamscan_annotation.out.kofamscan_tsv.collect()

        // combine kofamscan results
        combine_kofamscan_results(all_kofamscan_tsvs)
    }
}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '2 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv
    path("*.txt"), emit: stb_log, optional: true

    script:
    """
    python ${baseDir}/bin/generate_genome_stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    publishDir "${params.outdir}/predicted_orfs", mode: 'copy'


    errorStrategy 'ignore'
    
    memory = "6 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/pyrodigal:3.4.1--py310h4b81fae_0"

    input:
    tuple val(genome_name), path(fasta)

    output:
    tuple val(genome_name), path("*.ffn"), emit: predicted_orfs_ffn
    tuple val(genome_name), path("*.faa"), emit: predicted_orfs_faa
    tuple val(genome_name), path("*.gbk"), emit: predicted_orfs_gbk

    script:
    """
    pyrodigal \\
        -i ${fasta} \\
        -f "gbk" \\
        -o "${genome_name}.gbk" \\
        -d ${genome_name}.ffn \\
        -a ${genome_name}.faa \\
        --no-stop-codon
    """
}

process convert_gbk_to_gff {
    tag "${genome_name}_convert_gbk_to_gff"
    publishDir "${params.outdir}/predicted_orfs", mode: 'copy'


    errorStrategy 'ignore'
    
    memory = "2 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    tuple val(genome_name), path(gbk_file)

    output:
    tuple val(genome_name), path("*.gff"), emit: predicted_orfs_gff

    script:
    """
    python ${baseDir}/bin/gbk_to_gff.py --gbk ${gbk_file} --output ${genome_name}.gff
    """
}

process smorfinder_pre_called {
    tag "${genome_name}_smorfinder_pre_called"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    errorStrategy 'ignore'
    // rarely some genomes will fail for no discernible reason, skip over these

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:modf"

    input:
    tuple val(genome_name), path(fasta), path(gff), path(ffn), path(faa)

    output:
    path("*_smorfinder.gff"), emit: smorf_gff
    path("*_smorfinder.faa"), emit: smorf_faa
    path("*_smorfinder.ffn"), emit: smorf_ffn
    path("*_smorfinder.tsv"), emit: smorf_tsv

    script:
    """
    smorf pre-called ${fasta} ${faa} ${ffn} ${gff} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff ${genome_name}_smorfinder.gff
    ln -s ${genome_name}/${genome_name}.faa ${genome_name}_smorfinder.faa
    ln -s ${genome_name}/${genome_name}.ffn ${genome_name}_smorfinder.ffn
    ln -s ${genome_name}/${genome_name}.tsv ${genome_name}_smorfinder.tsv
    """
}

process smorfinder_single {
    tag "${genome_name}_smorfinder_single"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    errorStrategy 'ignore'
    // rarely some genomes will fail for no discernible reason, skip over these

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:modf"

    input:
    tuple val(genome_name), path(fasta)

    output:
    path("*.gff"), emit: smorf_gff
    path("*.faa"), emit: smorf_faa
    path("*.ffn"), emit: smorf_ffn
    path("*_smorfinder.tsv"), emit: smorf_tsv

    script:
    """
    smorf single ${fasta} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff
    ln -s ${genome_name}/${genome_name}.faa
    ln -s ${genome_name}/${genome_name}.ffn
    ln -s ${genome_name}/${genome_name}.tsv ${genome_name}_smorfinder.tsv
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

process filter_small_proteins {
    tag "${genome_name}_filter_small_proteins"

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"

    input:
    tuple val(genome_name), path(predicted_orfs_faa)

    output:
    tuple val(genome_name), path("*.fasta"), emit: filtered_proteins

    script:
    """
    python ${baseDir}/bin/filter_proteins.py ${predicted_orfs_faa} ${genome_name}_filtered_proteins.fasta --max_length 100
    """
}

process predict_cleavage_peptides {
    tag "${genome_name}_predict_cleavage_peptides"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    errorStrategy 'ignore'

    memory = '30 GB'
    cpus = 12

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/deeppeptide:v0.6.data"

    input:
    tuple val(genome_name), path(predicted_orfs_faa)

    output:
    tuple val(genome_name), path("${genome_name}_peptide_predictions.json"), emit: cleavage_peptides_json

    script:
    """
    WORKDIR=\$PWD
    cp ${predicted_orfs_faa} /app/DeepPeptide/predictor
    cd /app/DeepPeptide/predictor

    python3 predict.py --fastafile ${predicted_orfs_faa.getName()} --output_dir ${genome_name} --output_fmt json --batch_size 200
    
    cp ${genome_name}/peptide_predictions.json \$WORKDIR/${genome_name}_peptide_predictions.json
    """
}

process extract_cleavage_peptides_json {
    tag "${genome_name}_extract_cleavage_peptides_json"
    publishDir "${params.outdir}/cleavage_peptides", mode: 'copy'

    errorStrategy 'ignore'

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

    errorStrategy 'ignore'

    memory = "15 GB"
    cpus = 4

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

    memory = "20 GB"
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
    --output-counts all_molecule_counts.tsv \\
    --output-smorfinder all_smorfinder_results.tsv \\
    --output-deeppeptide all_deeppeptide_results.tsv
    """
}

process kofamscan_annotation {
    tag "${genome_name}_kofam_scan_annotation"
    publishDir "${params.outdir}/kofam_scan_annotation", mode: 'copy'

    errorStrategy 'ignore'

    memory = "15 GB"
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

    container "quay.io/biocontainers/polars:0.18.15"

    input:
    path(kofamscan_tsvs)

    output:
    path("combined_kofamscan_results.tsv"), emit: combined_kofamscan_results_tsv

    script:
    """
    python3 ${baseDir}/bin/combine_kofamscan_tsvs.py \\
    --input_files ${kofamscan_tsvs.join(' ')} \\
    --output combined_kofamscan_results.tsv
    """
}

