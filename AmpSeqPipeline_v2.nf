#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fastqDir  = 'fastq/'
params.outDir    = 'Results/'
params.threads   = 4
params.reference = params.reference ?: { error "Please provide --reference /path/to/genome.fa" }()
params.bed       = 'bed.txt'
params.tag       = 'Experiment01'

process FastQC {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/00_fastQC", mode: 'copy'
    cpus params.threads

    input:
    path fileList

    output:
    path ("*.html"), emit: html

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${fileList}
    """
}

process TrimNGSAdapters {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/01_trimmed", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("${sampleId}_R1_trimmed.fastq.gz"), path("${sampleId}_R2_trimmed.fastq.gz"), emit: trimmed
    path "${sampleId}_R{1,2}_unpaired.fastq.gz", emit: unpaired

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        ${r1} ${r2} \\
        ${sampleId}_R1_trimmed.fastq.gz ${sampleId}_R1_unpaired.fastq.gz \\
        ${sampleId}_R2_trimmed.fastq.gz ${sampleId}_R2_unpaired.fastq.gz \\
        ILLUMINACLIP:/usr/local/share/trimmomatic/TruSeq3-PE.fa:2:30:10:2:True \\
        LEADING:3
    """
}

process mergeReads {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/02_merged", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(r1), path(r2)

    output:
    tuple val(sampleId), path("${sampleId}_merged.fastq.gz"), emit: merged

    script:
    """
    flash ${r1} ${r2} \\
        --threads ${task.cpus} \\
        --max-overlap=320 \\
        --min-overlap=10 \\
        --output-prefix=${sampleId}_merged \\
        --compress

    mv "${sampleId}_merged.extendedFrags.fastq.gz" "${sampleId}_merged.fastq.gz"
    """
}

process mapBWA {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/03_mapped", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(mergedRead)

    output:
    tuple val(sampleId), path("${sampleId}.sorted.bam"), emit: mapped

    script:
    """
    bwa mem -t ${task.cpus} ${params.reference} ${mergedRead} \\
        | samtools view -bS \\
        | samtools sort -o "${sampleId}.sorted.bam"
    """
}

process indexBAM {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/03_mapped", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(mappedRead)

    output:
    tuple val(sampleId), path(mappedRead), path("${sampleId}.sorted.bam.bai"), emit: indexed

    script:
    """
    samtools index ${mappedRead}
    """
}

process countIndels {
    container 'recouto23/tsai_indels:1.0.0'

    publishDir "${params.outDir}/04_indels", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(mappedRead), path(indexedRead)

    output:
    tuple val(sampleId), path("${sampleId}_countIndelsReport.tsv")

    script:
    """
    python3 "$workflow.projectDir/Scripts/count_indels_integrations.py" \\
        --bed $workflow.projectDir/${params.bed} \\
        --ref ${params.reference} \\
        --bam ${mappedRead} \\
        --out ${sampleId}_countIndelsReport.tsv
    """
}

process generateIndelReport {
    container 'recouto23/ampseq:1.0.0'

    publishDir "${params.outDir}/04_indels", mode: 'copy'
    cpus params.threads

    input:
    path(fileList)

    output:
    path("${params.tag}_CombinedIndelReport.xlsx"), emit: report
    path("${params.tag}_ReadPlot.png"),              emit: readDepth

    script:
    """
    #!/usr/bin/env python3
    import pandas
    import glob
    import matplotlib.pyplot as plt
    import seaborn

    files = glob.glob("*countIndelsReport.tsv")
    dfs = [pandas.read_csv(f, sep='\\t', header=None) for f in files]
    combinedDF = pandas.concat(dfs, axis=0, ignore_index=True)
    combinedDF = combinedDF.rename({0:"File", 1:"Site", 2:"Chrom", 3:"Start", 4:"Stop",
                                    5:"Seq", 6:"Indels", 7:"Integrations", 8:"Total Reads"}, axis=1)
    combinedDF['Indel %'] = float(0)
    for a in range(len(combinedDF)):
        combinedDF.loc[a, 'Indel %'] = 100*combinedDF.loc[a, 'Indels']/combinedDF.loc[a, 'Total Reads']

    seaborn.stripplot(combinedDF, x="Site", y="Total Reads", hue="File", palette='BuPu')
    plt.legend(loc='upper right', bbox_to_anchor=[1.30, 1.00])
    plt.savefig("${params.tag}_ReadPlot.png")

    combinedDF.to_excel("${params.tag}_CombinedIndelReport.xlsx", index=False)
    """
}

process indelPlot {
    container 'recouto23/indel_plot:1.0.0'

    publishDir "${params.outDir}/04_indels/00_indelPlots", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sampleId), path(bam), path(bai)

    output:
    path("${sampleId}_indelPlot.png"), emit: indelPlot

    script:
    """
    python3 "$workflow.projectDir/Scripts/IndelPlot.py" \\
        --chrom "chr3" \\
        --cutSite 46373214 \\
        --buffer 50 \\
        --inputDirectory "$workflow.projectDir/${params.outDir}/03_mapped/*.bam" \\
        --tag "${sampleId}"
    """
}

workflow {
    fastqs = Channel.fromPath("${params.fastqDir}/*.fastq.gz")
    FastQC(fastqs)

    readCh = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz", checkIfExists: true)
    TrimNGSAdapters(readCh)

    mergeReads(TrimNGSAdapters.out.trimmed)
    mapBWA(mergeReads.out.merged)
    indexBAM(mapBWA.out.mapped)
    countIndels(indexBAM.out.indexed)

    generateIndelReport(countIndels.out.map { sampleId, tsv -> tsv }.collect())
    indelPlot(indexBAM.out.indexed)
}
