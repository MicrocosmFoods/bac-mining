#! /usr/bin/env nextflow

// Description
// Mine bacterial MAGs from fermented foods for peptides and BGCs

nextflow.enable.dsl=2

params.threads=16
params.outdir=null

log.info """\

MINE FERMENTED FOOD BACTERIAL GENOMES FOR BIOACTIVE PEPTIDES AND BGCs.

NOTE: YOU MUST PRE-DOWNLOAD THE ANTISMASH AND PFAM DATABASES AND PROVIDE THE PATH
WITH --antismash_db AND --pfam_db. THIS WORKFLOW DOES NOT SUPPORT DOWNLOADING 
DATABASES AUTOMATICALLY.
=================================================================
input_genomes                   : $params.input_genomes
genome_metadata                 : $params.genome_metadata
antismash_db                    : $params.antismash_db
pfam_db                         : $params.pfam_db
peptides_fasta                  : $params.peptides_fasta
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels and workflow steps
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [file, baseName]
    }

genome_metadata = channel.fromPath(params.genome_metadata)
antismash_db_ch = channel.fromPath(params.antismash_db)
peptides_db_ch = channel.fromPath(params.peptides_fasta)
pfam_db_ch = channel.fromPath(params.pfam_db)

workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[1] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    
    // get small ORF predictions with smorfinder
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file

    // combine smorf proteins into a single FASTA
    combine_smorf_proteins(smorf_proteins.collect())
    combined_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // cluster smorf proteins 100% identity and get representative seqs
    mmseqs_100id_cluster(combined_smorf_proteins)
    nonredundant_smorfs = mmseqs_95id_cluster.out.nonredundant_seqs_fasta
    mmseqs_clusters = mmseqs_95id_cluster.out.cluster_summary_tsv

    // mmseqs cluster summaries and stats
    summarize_mmseqs_clusters(mmseqs_clusters, nonredundant_smorfs, genome_metadata)

    // predict ORFs with pyrdogial
    pyrodigal(genome_fastas)
    predicted_orfs = pyrodigal.out.predicted_orfs_gbk

    // antismash predictions and extract info from GBKs
    antismash_input_ch = predicted_orfs.combine(antismash_db_ch)
    antismash(antismash_input_ch)
    antismash_gbk_files = antismash.out.gbk_results

    // bigscape on all antismash gbk_files
    all_antismash_gbk_files = antismash_gbk_files.map{ it[1] }.collect()
    run_bigscape(all_antismash_gbk_files, pfam_db_ch)
    bigscape_annotations_tsv = run_bigscape.out.bigscape_annotations_tsv

    // combine bigscape aggregate TSV with metadata
    combine_bigscape_metadata(bigscape_annotations_tsv, genome_metadata, genome_stb)

    // deepsig predictions on combined, non-redundant smorf proteins
    deepsig(nonredundant_smorfs)

    // peptides.py sequence characterization on combined, non-redundant smorf proteins
    characterize_peptides(nonredundant_smorfs)

    // DIAMOND seq similarity to Peptipedia peptide sequences of interest
    make_diamond_db(peptides_db_ch)
    peptides_dmnd_db = make_diamond_db.out.peptides_diamond_db
    diamond_blastp(nonredundant_smorfs, peptides_dmnd_db)

}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv

    script:
    """
    python ${baseDir}/bin/generate-genome-stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    memory = '10 GB'
    cpus = 4

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:latest"
    conda "envs/smorfinder.yml"

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
    ln -s ${genome_name}/${genome_name}.tsv
    """

}

process combine_smorf_proteins {
    tag "combine_smorf_proteins"
    publishDir "${params.outdir}/combined_smorf_proteins", mode: 'copy'

    memory = '10 GB'
    cpus = 1
    
    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(smorf_proteins)

    output: 
    path("*.fasta"), emit: combined_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} combined_smorf_proteins.fasta
    """
    
}

process mmseqs_95id_cluster {
    tag "mmseqs_95id_cluster"
    publishDir "${params.outdir}/mmseqs_95id_cluster", mode: 'copy'

    memory = '10 GB'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2"
    conda "envs/mmseqs2.yml"

    input:
    path(protein_fasta_file)
    
    output:
    path("*_rep_seq.fasta"), emit: nonredundant_seqs_fasta
    path("*.tsv"), emit: cluster_summary_tsv

    script:
    """
    mmseqs easy-cluster ${protein_fasta_file} nonredundant_smorf_proteins tmp --min-seq-id .95 --threads ${task.cpus}
    """   
}

process summarize_mmseqs_clusters {
    tag "summarize_mmseqs_clusters"
    publishDir "${params.outdir}/main_results/mmseqs_clusters", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(mmseqs_cluster_file)
    path(mmseqs_nonredundant_seqs)
    path(genome_metadata_tsv)

    output:
    path("mmseqs_summary.tsv"), emit: mmseqs_summary
    path("mmseqs_metadata.tsv"), emit: mmseqs_metadata
    path("mmseqs_substrate_counts.tsv"), emit: mmseqs_substrate_counts
    path("mmseqs_phylo_groups_counts.tsv"), emit: mmseqs_phylo_groups_counts

    script:
    """
    python ${baseDir}/bin/process_mmseqs_clusters.py ${mmseqs_cluster_file} ${mmseqs_nonredundant_seqs} ${genome_metadata_tsv} mmseqs_summary.tsv mmseqs_metadata.tsv mmseqs_substrate_counts.tsv mmseqs_phylo_groups_counts.tsv
    """ 
}

process pyrodigal {
    tag "${genome_name}_pyrodigal"
    
    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/pyrodigal:3.4.1--py310h4b81fae_0"
    conda "envs/pyrodigal.yml"

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
        -a ${genome_name}.faa
    """
}

process antismash {
    tag "${genome_name}_antismash"
    publishDir "${params.outdir}/antismash", mode: 'copy'

    memory = "20 GB"
    cpus = 4

    container "public.ecr.aws/biocontainers/antismash-lite:7.1.0--pyhdfd78af_0"
    conda "envs/antismashlite.yml"

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

process run_bigscape {
    tag "bigscape_all_gbks"
    publishDir "${params.outdir}/bigscape", mode: 'copy'

    memory = "20 GB"
    cpus = 6
    
    container "quay.io/biocontainers/bigscape:1.1.9--pyhdfd78af_0"
    conda "envs/bigscape.yml"

    input:
    path(gbk_files)
    path(pfam_db)

    output:
    path("*"), emit: bigscape_results
    path("network_files/*/Network_Annotations_Full.tsv"), emit: bigscape_annotations_tsv

    script:
    """
    bigscape -i ./ -o bigscape_results --pfam_dir ${pfam_db} --cores ${task.cpus} --mibig
    """

}

process combine_bigscape_metadata {
    tag "combine_bigscape_metadata"
    publishDir "${params.outdir}/main_results/bgc_info", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/csgenetics/tidyverse:latest"
    conda "envs/tidyverse.yml"

    input:
    path(bigscape_annotations_tsv)
    path(genome_metadata)
    path(genome_stb)

    output:
    path("bgc_annotation_metadata.tsv"), emit: bgc_metadata_tsv
    path("bgc_substrate_type_counts.tsv"), emit: bgc_substrates_tsv
    path("bgc_phylo_groups_counts.tsv"), emit: bgc_phylo_groups_tsv
    
    script:
    """
    Rscript ${baseDir}/bin/combine-bgc-metadata.R ${genome_metadata} ${bigscape_annotations_tsv} ${genome_stb} bgc_annotation_metadata.tsv bgc_substrate_type_counts.tsv bgc_phylo_groups_counts.tsv
    """
}

process deepsig {
    tag "deepsig_predictions"
    publishDir "${params.outdir}/deepsig", mode: 'copy'
    
    accelerator 1, type: 'nvidia-t4'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/deepsig:1.2.5--pyhca03a8a_1"
    conda "envs/deepsig.yml"

    input: 
    path(faa_file)

    output: 
    path("*.tsv"), emit: deepsig_tsv

    script: 
    """
    deepsig -f ${faa_file} -o nonredundant_smorf_proteins_deepsig.tsv -k gramp -t ${task.cpus}
    """

}

process characterize_peptides {
    tag "characterize_peptides"
    publishDir "${params.outdir}/peptide_characterization", mode: 'copy'

    memory = "5 GB"
    cpus = 1

    container "elizabethmcd/peptides"
    conda "envs/peptides.yml"

    input:
    path(faa_file)

    output: 
    path("*.tsv"), emit: peptides_tsv

    script:
    """
    python ${baseDir}/bin/characterize_peptides.py ${faa_file} nonredundant_smorf_proteins_peptide_characteristics.tsv
    """
}

process make_diamond_db {
    tag "make_diamond_db"

    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/diamond:2.1.7--h43eeafb_1"
    conda "envs/diamond.yml"

    input:
    path(peptides_fasta)

    output:
    path("*.dmnd"), emit: peptides_diamond_db

    script:
    """
    diamond makedb --in ${peptides_fasta} -d peptides_db.dmnd
    """
}

process diamond_blastp {
    tag "diamond_blastp"
    publishDir "${params.outdir}/diamond_blastp", mode: 'copy'

    memory = "10 GB"

    container "public.ecr.aws/biocontainers/diamond:2.1.7--h43eeafb_1"
    conda "envs/diamond.yml"

    input:
    path(faa_file)
    path(peptides_diamond_db)

    output:
    path("*.tsv"), emit: blastp_hits_tsv

    script:
    """
    diamond blastp -d ${peptides_diamond_db} -q ${faa_file} -o nonredundant_smorf_proteins_blast_results.tsv --header simple \\
    --outfmt 6 qseqid sseqid full_sseq pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
    """
}
