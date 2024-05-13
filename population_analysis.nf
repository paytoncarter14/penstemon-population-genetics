params.reference_fasta = "/nobackup/scratch/grp/fslg_pws670/payton/p_eatonii_hard_masked.fasta"
params.bams_dir = "${projectDir}/bams"
params.populations_dir = "${projectDir}/populations"
params.comparisons = [['eatonii', 'exsertus'], ['eatonii', 'undosus'], ['exsertus', 'undosus'], ['Utah', 'Arizona'], ['all']]

// file prep

process get_bam_paths {
    input: path(sample_list)
    output: tuple path(sample_list), path("*.bamPaths.txt")

    executor 'local'

    shell: '''sed -e 's#^#!{params.bams_dir}/#g' -e 's#$#.bam#g' !{sample_list} > !{sample_list.simpleName}.bamPaths.txt'''
    stub: "touch ${sample_list.simpleName}.bamPaths.txt"
}

// saf, sfs, and thetas

process saf {
    input: tuple path(sample_list), path(bam_paths)
    output: tuple val(bam_paths.simpleName), path("*.saf.idx"), path("*.saf.pos.gz"), path("*.saf.gz")

    cpus 8
    memory {bam_paths.simpleName == "all" ? 300.GB + (100.GB * task.attempt) : 100.GB * task.attempt}
    time {bam_paths.simpleName == "all" ? 2.hour * task.attempt : 8.hour * task.attempt}

    shell: '''angsd -b !{bam_paths} -anc !{params.reference_fasta} -out !{bam_paths.simpleName} -doSaf 1 -GL 1 -nThreads !{task.cpus}'''
    stub: "x=${sample_list.simpleName}; touch \${x}.saf.idx \${x}.saf.pos.gz \${x}.saf.gz"
}

process sfs {
    input: tuple val(population), path(saf_idx), path(saf_pos_gz), path(saf_gz)
    output: tuple val(population), path("*.sfs")

    cpus 8
    memory {32.GB * task.attempt}
    time {1.hour * task.attempt}

    shell: '''realSFS !{saf_idx} -fold 1 -P !{task.cpus} > !{population}.sfs'''
    stub: "touch ${population}.sfs"
}

process thetas {
    input: tuple val(population), path(saf_idx), path(saf_pos_gz), path(saf_gz), path(sfs)
    output: tuple val(population), path("*.thetas.gz"), path("*.thetas.idx")

    cpus 8
    memory {8.GB * task.attempt}
    time {1.hour * task.attempt}

    shell: '''realSFS saf2theta !{saf_idx} -sfs !{sfs} -outname !{population} -P !{task.cpus}'''
    stub: "touch ${population}.thetas.gz ${population}.thetas.idx"
}

process thetaStat {
    input: tuple val(population), path(thetas_gz), path(thetas_idx)
    output: tuple val(population), path("*.pestPG")

    cpus 1
    memory {8.GB * task.attempt}
    time {2.hour * task.attempt}

    shell: '''thetaStat do_stat !{thetas_idx} -win 50000 -step 10000 -outnames !{population}'''
    stub: "touch ${population}.pestPG"
}

// pca and admixture

process combine_sample_lists {
    input: val(comparison)
    output: tuple path("*.bamPaths.txt"), path("*.key.txt")

    executor 'local'

    shell:
    if (comparison.size == 1) {
        '''
        sed -e 's#^#!{params.bams_dir}/#g' -e 's#$#.bam#g' !{params.populations_dir}/!{comparison[0]}.txt > !{comparison[0]}.bamPaths.txt
        sed -e 's#$# !{comparison[0]}#g' !{params.populations_dir}/!{comparison[0]}.txt > !{comparison[0]}.key.txt
        '''
    } else {
        fn = "${comparison[0]}_v_${comparison[1]}"
        '''
        cat !{params.populations_dir}/!{comparison[0]}.txt !{params.populations_dir}/!{comparison[1]}.txt | sed -e 's#^#!{params.bams_dir}/#g' -e 's#$#.bam#g' > !{fn}.bamPaths.txt
        sed -e 's#$# !{comparison[0]}#g' !{params.populations_dir}/!{comparison[0]}.txt >> !{fn}.key.txt
        sed -e 's#$# !{comparison[1]}#g' !{params.populations_dir}/!{comparison[1]}.txt >> !{fn}.key.txt
        '''
    }
}

process gl {
    input: tuple path(bam_paths), path(key)
    output: tuple path(key), path("*.beagle.gz"), path("*.mafs.gz")

    cpus 8
    memory {180.GB * task.attempt}
    time {3.hour * task.attempt}

    shell: '''angsd -GL 1 -out !{key.simpleName} -nThreads !{task.cpus} -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam !{bam_paths}'''
    stub: "x=${key.simpleName}; touch \${x}.beagle.gz \${x}.mafs.gz"
}

process pcangsd {
    input: tuple path(key), path(beagle_gz), path(mafs_gz)
    output: tuple path(key), path("*.cov")

    cpus 8
    memory 32.GB
    time 1.hour

    shell: '''python3 !{projectDir}/bin/pcangsd-v.0.99/pcangsd.py -beagle !{beagle_gz} -o !{beagle_gz.simpleName}_pca -threads !{task.cpus}'''
    stub: "touch ${beagle_gz.simpleName}.cov"
}

process ngsadmix {
    input:
        tuple path(key), path(beagle_gz), path(mafs_gz)
        each k
    output: tuple path(key), path("*.qopt"), path("*.fopt.gz"), path("*.log")

    cpus 8
    memory {32.GB * task.attempt}
    time {4.hour * task.attempt}

    shell: '''NGSadmix -likes !{beagle_gz} -K !{k} -o !{beagle_gz.simpleName}_admix_k!{k} -P !{task.cpus}'''
    stub: "x=${beagle_gz.simpleName}; touch \${x}.qopt \${x}.fopt.gz \${x}.log"
}

// psmc

process mpileup {
    input: path(bam)
    output: path("*.diploid.fq.gz")

    cpus 1
    memory {8.GB * task.attempt}
    time {3.hour * task.attempt}

    shell: '''
    sample=$(basename !{bam} .bam)
    coverage=$(grep -P "${sample}\\t" !{params.bams_dir}/coverage.txt | cut -f 2)
    min_cov=$(echo "${coverage}/3" | bc)
    max_cov=$(echo "${coverage}*2" | bc)
    echo "hi"

    samtools mpileup -C50 -uf !{params.reference_fasta} !{bam} | bcftools call -c - | vcfutils.pl vcf2fq -d ${min_cov} -D ${max_cov} | gzip > !{bam.simpleName}.diploid.fq.gz
    '''
    stub: "touch ${bam.simpleName}.diploid.fq.gz"
}

process fq2psmcfa {
    input: path(diploid_fq_gz)
    output: path("*.psmcfa")

    cpus 1
    memory 4.GB
    time 15.minute

    shell: '''fq2psmcfa -q20 !{diploid_fq_gz} > !{diploid_fq_gz.simpleName}.psmcfa'''
    stub: "touch ${diploid_fq_gz.simpleName}.psmcfa"
}

process psmc {
    input: path(psmcfa)
    output: path("*.psmc")

    cpus 1
    memory 4.GB
    time 15.minute

    shell: '''psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o !{psmcfa.simpleName}.psmc !{psmcfa}'''
    stub: "touch ${psmcfa.simpleName}.psmc"
}

process psmc_plot {
    // peatonii594_diploid.psmc gives a divide by zero error when included, lots of NaN in psmc file
    input: path(psmc_files)
    output: path("*.pdf")

    executor 'local'

    shell: '''psmc_plot.pl -p -u 7.0e-09 -g 1 psmc !{psmc_files}'''
    stub: "touch psmc.pdf"
}

workflow {
    // saf, sfs, and thetas
    get_bam_paths_c = get_bam_paths(Channel.fromPath("${params.populations_dir}/*.txt"))
    saf_c = saf(get_bam_paths_c)
    sfs_c = sfs(saf_c)
    thetas_c = thetas(saf_c.join(sfs_c))
    thetaStat_c = thetaStat(thetas_c)

    // pca and admixture
    combine_sample_lists_c = combine_sample_lists(Channel.fromList(params.comparisons))
    gl_c = gl(combine_sample_lists_c)
    pcangsd_c = pcangsd(gl_c)
    ngsadmix_c = ngsadmix(gl_c, [2, 3, 4])

    // psmc
    mpileup_c = mpileup(Channel.fromPath("${params.bams_dir}/*.bam"))
    fq2psmcfa_c = fq2psmcfa(mpileup_c)
    psmc_c = psmc(fq2psmcfa_c)
    psmc_plot_c = psmc_plot(psmc_c.collect())
}