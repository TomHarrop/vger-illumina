#!/usr/bin/env python3

import pathlib


#############
# FUNCTIONS #
#############

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
        expand(('output/04_meraculous/{read_set}_k{k}_diplo{diplo}/'
                'meraculous_final_results/final.scaffolds.fa'),
               read_set=['norm', 'trim-decon'],
               k=['31', '71', '127'],
               diplo=['0', '1'])

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
        dmin = '21'
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
            'echo \''
            'bin/meraculous/run_meraculous.sh '
            '-dir {params.outdir} '
            '-config {output.config} '
            '&> {log}'
            '\'')

# 03 normalise input
rule plot_kmer_coverage:
    input:
        hist_before = 'output/030_norm/hist.txt',
        hist_after = 'output/030_norm/hist_out.txt',
        peaks = 'output/030_norm/peaks.txt'
    output:
        plot = 'output/030_norm/kmer_plot.pdf'
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
        norm = 'output/logs/030_norm.log'
    threads:
        50
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target=60 '
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
        20
    shell:
        'bbmerge.sh '
        'threads={threads} '
        'in={input.fq} '
        'verystrict=t '
        'out={output.merged} '
        'outu={output.unmerged} '
        'ihist={output.ihist} '
        '2> {log.merge} '

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
        20
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
