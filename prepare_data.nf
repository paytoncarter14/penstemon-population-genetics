/*

fastp 0.23.2
multiqc 1.19
minimap2 2.27-r1193
samtools 1.16.1

hard mask assembly
sed '/^>/!s/[atgc]/N/g' p-eatonii_1252744.review.fasta.renamed.masked.no-cp-scaffolds.fasta > p_eatonii_hard_masked.fasta

to get list of samples:
find /nobackup/scratch/grp/fslg_pws670/penstemon_data/Illumina_reseq \
/nobackup/scratch/grp/fslg_pws670/pws672_w2024/2_resequencing_mapping/reseq_carrie \
-iname "*.fastq.gz" -not -iname "*trimmed_reads*" -exec readlink -e {} \; | grep -v "trimmed_reads\|test_r"

created symlinks in fastq directory and manually renamed links so all had [sample_id]_L[]_R[]_001.fastq.gz

looks like 26 samples? should be 20 from us, 6 from Carrie Wessinger lab. 36 total files. some samples have multiple sequencing runs and will need to be concatenated
ls /nobackup/scratch/grp/fslg_pws670/payton/fastq | cut -d '_' -f 1 | sort | uniq | wc -l

added variety and GPS coords to sampling spreadsheet

indexed hard masked assembly
bin/minimap2 -x sr -d p_eatonii_hard_masked.mmi p_eatonii_hard_masked.fasta

*/

run_name = 'unmasked'

index_path = '/nobackup/scratch/grp/fslg_pws670/payton/penstemon_final_nuclear_cp_assembly.mmi'
assembly_path = '/nobackup/scratch/grp/fslg_pws670/payton/penstemon_final_nuclear_cp_assembly.fasta'
fastq_dir = '/nobackup/scratch/grp/fslg_pws670/payton/all_fastq'

out_dir = 'output/' + run_name

process fastp {
    input: tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("${sample_id}_trimmed_R?.fastq.gz"), path("${sample_id}_fastp.json"), path("${sample_id}_fastp.html")
    publishDir path: out_dir + '/fastp/fastq', pattern: '*.fastq.gz'
    publishDir path: out_dir + '/fastp/html', pattern: '*.html'
    publishDir path: out_dir + '/fastp/json', pattern: '*.json'
    time '2h'
    memory '8 GB'
    shell: '''
    fastp \
    -i !{reads[0]} \
    -I !{reads[1]} \
    -o !{sample_id}_trimmed_R1.fastq.gz \
    -O !{sample_id}_trimmed_R2.fastq.gz \
    -j !{sample_id}_fastp.json \
    -h !{sample_id}_fastp.html \
    -Q -g
    '''
}

process multiqc {
    input: path '*_fastp.json'
    output: path 'multiqc_report.html'
    publishDir out_dir + '/fastp'
    time '10 m'
    memory '4 GB'
    shell: 'multiqc .'
}

process minimap2 {
    input: tuple val(sample_id), path(reads)
    output: tuple path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    publishDir out_dir + '/minimap2', mode: 'copy', overwrite: 'true'
    time '24 h'
    memory '16 GB'
    cpus '3'
    shell: '''
    minimap2 -ax sr !{index_path} !{reads[0]} !{reads[1]} | samtools sort -o !{sample_id}.bam - 
    samtools index !{sample_id}.bam
    '''
}

process flagstat {
    input: tuple path(bam), path(bai)
    output: path '*.bam.flagstat'
    publishDir out_dir + '/minimap2'
    time '1 h'
    memory '8 GB'
    shell: 'samtools flagstat !{bam} > !{bam}.flagstat'
}

process blobtools {
    input: tuple path(bam), path(bai)
    output: path '*.bam.cov'
    publishDir out_dir + '/minimap2'
    time '1 h'
    memory '8 GB'
    shell: 'blobtools map2cov -i !{assembly_path} -b !{bam}'
}

process get_stats {
    input: path '*'
    output: path 'stats.csv'
    publishDir out_dir
    time '15 m'
    memory '8 GB'
    shell:
    '''
    get_stats.py
    '''
}

workflow {

    init_c = Channel.fromFilePairs(fastq_dir + '/*_R{1,2}.fastq.gz')
    fastp_c = fastp(init_c)
    multiqc_c = multiqc(fastp_c.map{it[2]})
    minimap2_c = minimap2(fastp_c.map{it[0..1]})
    flagstat_c = flagstat(minimap2_c)
    blobtools_c = blobtools(minimap2_c)
    get_stats_c = get_stats(fastp_c.map{it[2]}.concat(flagstat_c).concat(blobtools_c).collect())

}