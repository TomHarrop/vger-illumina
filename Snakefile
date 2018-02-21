#!/usr/bin/env python3


###########
# GLOBALS #
###########

r1_raw = 'data/PM-OS170002-07/rawdata/VB_R1.fq.gz'
r2_raw = 'data/PM-OS170002-07/rawdata/VB_R2.fq.gz'
bbduk_ref = 'venv/bin/resources/phix174_ill.ref.fa.gz'
bbduk_adaptors = 'venv/bin/resources/adapters.fa'

#########
# SETUP #
#########

#########
# RULES #
#########

rule target:
    input:
        'output/010_trim-decon/vger.fq.gz'

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
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'bin/bbmap/repair.sh '
        'in=stdin.fastq ' 
        'out={output.fq} '
        '2> {log.repair} '
