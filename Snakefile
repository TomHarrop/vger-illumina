#!/usr/bin/env python3

import os
import pathlib


#############
# FUNCTIONS #
#############

def find_completed_assemblies(wildcards):
    my_files = list((dirpath, filenames)
                    for (dirpath, dirnames, filenames) 
                    in os.walk('output/04_meraculous'))
    my_fasta_files = []
    for dirpath, filenames in my_files:
        for filename in filenames:
            if ('final.scaffolds.fa' in filename 
                    and 'meraculous_final_results' in dirpath):
                my_path = os.path.join(dirpath, filename)
                my_fasta_files.append(my_path)
    return(my_fasta_files)


def resolve_path(x):
    mypath = pathlib.Path(x).resolve()
    return str(mypath)


def readset_wildard_resolver(wildcards):
    if wildcards.read_set == 'norm':
        return({'fq': 'output/030_norm/vger.fq.gz'})
    elif wildcards.read_set == 'trim-decon':
        return({'fq': 'output/010_trim-decon/vger.fq.gz'})
    else:
        raise ValueError('unknown read_set')


###########
# GLOBALS #
###########

r1_raw = 'data/PM-OS170002-07/rawdata/VB_R1.fq.gz'
r2_raw = 'data/PM-OS170002-07/rawdata/VB_R2.fq.gz'
bbduk_ref = 'venv/bin/resources/phix174_ill.ref.fa.gz'
bbduk_adaptors = 'venv/bin/resources/adapters.fa'
meraculous_config_file = 'src/meraculous_config.txt'

# containers
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'


#########
# SETUP #
#########

# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())

#########
# RULES #
#########

rule target:
    input:
        'output/011_kraken/kraken_report.txt',
        'output/020_merge/ihist.txt',
        'output/030_norm/kmer_plot.pdf',
        expand(('output/04_meraculous/{read_set}_k{k}_diplo{diplo}/'
                'meraculous_final_results/final.scaffolds.fa'),
               read_set=['norm', 'trim-decon'],
               k=['31', '71', '101'],
               diplo=['0', '1']),
        expand(('output/04_meraculous/trim-decon_k{k}_diplo0/'
                'meraculous_final_results/final.scaffolds.fa'),
               k=['61', '67', '69', '73', '75', '81', '91'])


# 05 run bbmap stats on completed assemblies
rule stats_plot:
    input:
        stats = 'output/050_assembly-stats/stats.txt'
    output:
        plot = 'output/050_assembly-stats/assembly_stats.pdf'
    log:
        log = 'output/logs/050_assembly-stats/plot.log'
    script:
        'src/plot_assembly_stats.R'

rule bbmap_stats:
    input:
        fa = find_completed_assemblies
    output:
        stats = 'output/050_assembly-stats/stats.txt'
    log:
        'output/logs/050_assembly-stats/stats.log'
    threads:
        1
    run:
        my_inputfiles = ','.join(input.fa)
        shell('statswrapper.sh '
              'in={my_inputfiles} '
              'minscaf=1000 '
              'format=3 '
              '> {output.stats} '
              '2> {log}')

# 04 launch meraculous
rule meraculous:
    input:
        unpack(readset_wildard_resolver)
    output:
        config = ('output/04_meraculous/'
                  '{read_set}_k{k}_diplo{diplo}/config.txt'),
        contigs = ('output/04_meraculous/'
                   '{read_set}_k{k}_diplo{diplo}/'
                   'meraculous_final_results/final.scaffolds.fa'),
    params:
        outdir = 'output/04_meraculous/{read_set}_k{k}_diplo{diplo}/',
        dmin = '0'
    threads:
        50
    log:
        'output/logs/04_meraculous/{read_set}_k{k}_diplo{diplo}.log'
    run:
        my_fastq = resolve_path(input.fq)
        my_conf = meraculous_config_string.format(
            my_fastq,
            wildcards.k,
            wildcards.diplo,
            params.dmin,
            threads)
        with open(output.config, 'wt') as f:
            f.write(my_conf)
        shell(
            'run_meraculous.sh '
            '-dir {params.outdir} '
            '-config {output.config} '
            '-cleanup_level 2 '
            '&> {log}')

# 03 normalise input
rule plot_kmer_coverage:
    input:
        hist_before = 'output/030_norm/hist.txt',
        hist_after = 'output/030_norm/hist_out.txt',
        peaks = 'output/030_norm/peaks.txt'
    output:
        plot = 'output/030_norm/kmer_plot.pdf'
    threads:
        1
    log:
        log = 'output/logs/030_norm/plot_kmer_coverage.log'
    singularity:
        r_container
    script:
        'src/plot_kmer_coverage.R'

rule bbnorm:
    input:
        fq = 'output/010_trim-decon/vger.fq.gz'
    output:
        fq_norm = 'output/030_norm/vger.fq.gz',
        fq_toss = 'output/030_norm/toss.fq.gz',
        hist = 'output/030_norm/hist.txt',
        hist_out = 'output/030_norm/hist_out.txt',
        peaks = 'output/030_norm/peaks.txt'
    log:
        norm = 'output/logs/030_norm/bbnorm.log'
    params:
        target = 60
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min=5 '
        'peaks={output.peaks} '
        '2> {log.norm} '  

# 02 attempt to merge overlapping reads
rule bbmerge:
    input:
        fq = 'output/010_trim-decon/vger.fq.gz'
    output:
        merged = 'output/020_merge/merged.fq.gz',
        unmerged = 'output/020_merge/unmerged.fq.gz',
        ihist = 'output/020_merge/ihist.txt'
    log:
        merge = 'output/logs/020_merge.log'
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbmerge.sh '
        'threads={threads} '
        'in={input.fq} '
        'verystrict=t '
        'out={output.merged} '
        'outu={output.unmerged} '
        'ihist={output.ihist} '
        '2> {log.merge} '

# 011 run kraken on decontaminated reads
rule kraken:
    input:
        r1 = 'output/011_kraken/r1.fq.gz',
        r2 = 'output/011_kraken/r2.fq.gz',
        db = directory('data/20180917-krakendb')
    output:
        out = 'output/011_kraken/kraken_out.txt',
        report = 'output/011_kraken/kraken_report.txt'
    log:
        'output/logs/011_kraken/kraken.log'
    threads:
        50
    priority:
        10
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

rule split:
    input:
        fq = 'output/010_trim-decon/vger.fq.gz'
    output:
        r1 = temp('output/011_kraken/r1.fq.gz'),
        r2 = temp('output/011_kraken/r2.fq.gz')
    threads:
        25
    log:
        'output/logs/011_kraken/reformat.log'
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input.fq} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        '2> {log}'


# 01 trim and decontaminate reads
rule trim_decon:
    input:
        r1 = r1_raw,
        r2 = r2_raw
    output:
        fq = 'output/010_trim-decon/vger.fq.gz',
        f_stats = 'output/010_trim-decon/vger_filter-stats.txt',
        t_stats = 'output/010_trim-decon/vger_trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/filter.log',  
        trim = 'output/logs/010_trim-decon/trim.log',
        repair = 'output/logs/010_trim-decon/repair.log'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq ' 
        'out={output.fq} '
        '2> {log.repair} '
