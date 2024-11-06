#!/usr/bin/env snakemake

#original author: Shaghayegh Soudi
#contributors: NA


configfile:
    "config.yaml"

SAMPLES = glob_wildcards(config['data']+"/{sample}.fastq")

### demultiplexing
rule demultiplexing
    input:
        "indexes/{genome}/{genome}.fa"
    output:
        directory("indexes/{genome}/Bisulfite_Genome")
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        script
    wrapper:
        "v2.6.0/bio/bismark/bismark_genome_preparation"


### fastQC
rule run_fastQC_raw_fastq
    input:
        "indexes/{genome}/{genome}.fa"
    output:
        directory("indexes/{genome}/Bisulfite_Genome")
    log:
        "logs/indexes/{genome}/Bisulfite_Genome.log"
    params:
        ""  # optional params string
    wrapper:
        "v2.6.0/bio/bismark/bismark_genome_preparation"


rule trim_galore:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "trimmed_data/{sample}_R1_val_1.fq.gz",
        "trimmed_data/{sample}_R2_val_2.fq.gz",
        "trimmed_data/{sample}_trimming_report.txt"
    params:
        rrbs = "--rrbs",  # Trim Galore parameter for RRBS-specific trimming
    log:
        "logs/trim_galore/{sample}_trim.log"
    shell:
        """
        trim_galore {params.rrbs} --paired {params.adapter} \
        -o trimmed_data \
        {input[0]} {input[1]} > {log} 2>&1
        """


rule bismark_rrbs_alignment:
    input:
        fq_1="reads/{sample}.1.fastq",
        fq_2="reads/{sample}.2.fastq",
        genome="indexes/{genome}/{genome}.fa",
        bismark_indexes_dir="indexes/{genome}/Bisulfite_Genome",
        genomic_freq="indexes/{genome}/genomic_nucleotide_frequencies.txt"
    output:
        bam="bams/{sample}_{genome}_pe.bam",
        report="bams/{sample}_{genome}_PE_report.txt",
        nucleotide_stats="bams/{sample}_{genome}_pe.nucleotide_stats.txt",
        bam_unmapped_1="bams/{sample}_{genome}_unmapped_reads_1.fq.gz",
        bam_unmapped_2="bams/{sample}_{genome}_unmapped_reads_2.fq.gz",
        ambiguous_1="bams/{sample}_{genome}_ambiguous_reads_1.fq.gz",
        ambiguous_2="bams/{sample}_{genome}_ambiguous_reads_2.fq.gz"
    params:
        genome_dir="path/to/bismark_genome",  # Path to the Bismark-prepared genome directory
        extra="--rrbs --bowtie2 --unmapped"  # RRBS-specific, Bowtie2, and keep unmapped options
        params:
        # optional params string, e.g: -L32 -N0 -X400 --gzip
        # Useful options to tune:
        # (for bowtie2)
        # -N: The maximum number of mismatches permitted in the "seed", i.e. the first L base pairs
        # of the read (deafault: 1)
        # -L: The "seed length" (deafault: 28)
        # -I: The minimum insert size for valid paired-end alignments. ~ min fragment size filter (for
        # PE reads)
        # -X: The maximum insert size for valid paired-end alignments. ~ max fragment size filter (for
        # PE reads)
        # --gzip: Gzip intermediate fastq files
        # --ambiguous --unmapped
        # -p: bowtie2 parallel execution
        # --multicore: bismark parallel execution
        # --temp_dir: tmp dir for intermediate files instead of output directory
        extra=' --ambiguous --unmapped --nucleotide_coverage',
        basename='{sample}_{genome}'
    log:
        "logs/bismark/{sample}_alignment.log"
    threads: 8  # Adjust based on computational resources
    shell:
        """
        bismark {params.genome_dir} -1 {input.fastq1} -2 {input.fastq2} \
        {params.extra} -p {threads} -o aligned_data > {log} 2>&1
        """


rule deduplicate_bam_rrbs:
    input:
        bam="aligned_data/{sample}_bismark_bt2.bam",
        bai="aligned_data/{sample}_bismark_bt2.bam.bai"  # BAM index file
    output:
        dedup_bam="dedup_data/{sample}_deduplicated.bam",
        dedup_log="dedup_data/{sample}_deduplication_report.txt"
    params:
        extra="--rrbs"  # RRBS-specific deduplication option, if applicable
    log:
        "logs/deduplication/{sample}_deduplication.log"
    shell:
        """
        deduplicate_BAM_files {params.extra} --input {input.bam} \
        --output {output.dedup_bam} > {output.dedup_log} 2> {log}
        """


rule bismark_methylation_extractor:
    input: "bams/{sample}.bam"
    output:
        mbias_r1="qc/meth/{sample}.M-bias_R1.png",
        # Only for PE BAMS:
        # mbias_r2="qc/meth/{sample}.M-bias_R2.png",

        mbias_report="meth/{sample}.M-bias.txt",
        splitting_report="meth/{sample}_splitting_report.txt",

        # 1-based start, 1-based end ('inclusive') methylation info: % and counts
        methylome_CpG_cov="meth_cpg/{sample}.bismark.cov.gz",
        # BedGraph with methylation percentage: 0-based start, end exclusive
        methylome_CpG_mlevel_bedGraph="meth_cpg/{sample}.bedGraph.gz",

        # Primary output files: methylation status at each read cytosine position: (extremely large)
        read_base_meth_state_cpg="meth/CpG_context_{sample}.txt.gz",
        # * You could merge CHG, CHH using: --merge_non_CpG
        read_base_meth_state_chg="meth/CHG_context_{sample}.txt.gz",
        read_base_meth_state_chh="meth/CHH_context_{sample}.txt.gz"
    log:
        "logs/meth/{sample}.log"
    params:
        output_dir="meth",  # optional output dir
        extra="--gzip --comprehensive --bedGraph"  # optional params string
    wrapper:
        "v2.6.0/bio/bismark/bismark_methylation_extractor"