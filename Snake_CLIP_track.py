#configfile: "config.yaml"
# snakemake -j 30 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
# snakemake -j 30 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn=4 -q home-yeo" 

import pandas as pd
import os
import sys
import glob
from snake_config import CHROM_SIZES, GENOME_FA, PYTHON3_PATH


include: "snake_config.py"

if not os.path.exists(MANIFEST): make_meta(MANIFEST)
manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')

sample_labels = manifest.Sample.tolist()

rule all:
    input:
        expand("CITS/{sample_label}.plus.bw", sample_label = sample_labels)+expand("CIMS/{sample_label}.plus.bw", sample_label = sample_labels)+expand("coverage/{sample_label}.plus.bw", sample_label = sample_labels)
    output:
        "snakeCLIP.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "1",
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        "echo $(date)  > {output};"
        "echo created by Evan Boyle and the Yeo lab >> {output}"


# this might not work for single-end?
# yes when single end we only have header
rule extract_read_two:
    input:
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
    output:
        read2="processed_bam/{sample_label}.r2.bam",
        read1="processed_bam/{sample_label}.r1.bam"
    params:
        run_time=6,
        error_out_file = "error_files/extract_read2",
    shell:
        """
        module load samtools
        paired=$(samtools view -c -f 1 {input.bam})
        if [ "$paired" -lt "1" ]; 
            then cp {input.bam} {output.read2};
            cp {input.bam} {output.read1}
        else


            # # get read2 only
            samtools view -h -f 0x0080 {input.bam} | samtools view -Sb - > {output.read2}

            # get read1 only
            samtools view -h -f 0x0040 {input.bam} | samtools view -Sb - > {output.read1}
        fi
        """
rule strand_specific_bam:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        forward_bam="processed_bam/{sample_label}.r2.forward.bam",
        reverse_bam="processed_bam/{sample_label}.r2.reverse.bam"
    params:
        run_time=6,
        error_out_file = "error_files/pileup",
    shell:
        """
        # make strand specific
        module load samtools
        
        samtools view -h -F 0x10 {input.bam} | samtools view -Sb - > {output.forward_bam}
        samtools view -h -f 0x10 {input.bam}| samtools view -Sb - > {output.reverse_bam}
        """


rule bam_pileup:
    input:
        forward_bam="processed_bam/{sample_label}.r2.forward.bam",
        reverse_bam="processed_bam/{sample_label}.r2.reverse.bam"
    output:
        plus_pileup="CIMS/{sample_label}.plus.pileup",
        minus_pileup="CIMS/{sample_label}.minus.pileup",
        
    params:
        run_time=6,
        genomefa=GENOME_FA,
        error_out_file = "error_files/bam_pileup",
    shell:
        """
        module load samtools
        
        samtools mpileup -s  {input.forward_bam} -f {params.genomefa} > {output.plus_pileup}
        samtools mpileup -s  {input.reverse_bam} -f {params.genomefa} > {output.minus_pileup}
        module unload samtools
        """
rule CIMS_bw:
    input:
        plus_pileup="CIMS/{sample_label}.plus.pileup",
        minus_pileup="CIMS/{sample_label}.minus.pileup",
    output:
        plus="CIMS/{sample_label}.plus.bw",
        minus="CIMS/{sample_label}.minus.bw"
    params:
        run_time=6,
        chr_size=CHROM_SIZES,
        python=PYTHON3_PATH,
        error_out_file = "error_files/CIMS_bw",

    shell:
        """
        # coverage
        {params.python} pileupToMismatchBw.py {input.plus_pileup} {params.chr_size} {output.plus}
        {params.python} pileupToMismatchBw.py {input.minus_pileup} {params.chr_size} {output.minus}
        """
rule CITS_bam_to_bedgrah:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        plus="CITS/{sample_label}.plus.bedgraph",
        minus="CITS/{sample_label}.minus.bedgraph"
    params:
        run_time=6,
        error_out_file = "error_files/CITS_bedgraph"
    shell:
        """
        module load bedtools;
        set +o pipefail;
        
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + -5 > {output.plus}
        bedSort {output.plus} {output.plus}

        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand - -5 > {output.minus}
        bedSort {output.minus} {output.minus}
        """
rule CITS_bedgraph_to_bw:
    input:
        plus="CITS/{sample_label}.plus.bedgraph",
        minus="CITS/{sample_label}.minus.bedgraph"
    output:
        plus="CITS/{sample_label}.plus.bw",
        minus="CITS/{sample_label}.minus.bw"
    params:
        run_time=14,
        chr_size=CHROM_SIZES,
        error_out_file = "error_files/CITS_bw",
    shell:
        """
        module load ucsctools
        bedGraphToBigWig {input.plus} {params.chr_size} {output.plus}
        bedGraphToBigWig {input.minus} {params.chr_size} {output.minus}
        """
rule COV_bam_to_bigwig_coverage:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        plus="coverage/{sample_label}.plus.bedgraph",
        minus="coverage/{sample_label}.minus.bedgraph"
    params:
        run_time=6,
        error_out_file = "error_files/coverage_bedgraph",
    shell:
        """
        module load bedtools;
        set +o pipefail;

        # coverage
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + > {output.plus}
        bedSort {output.plus} {output.plus}
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand - > {output.minus}
        bedSort {output.minus} {output.minus}
        """
rule COV_bedgraph_to_bw:
    input:
        plus="coverage/{sample_label}.plus.bedgraph",
        minus="coverage/{sample_label}.minus.bedgraph"
    output:
        plus="coverage/{sample_label}.plus.bw",
        minus="coverage/{sample_label}.minus.bw"
    params:
        run_time=14,
        chr_size=CHROM_SIZES,
        error_out_file = "error_files/CITS_bw",
    shell:
        """
        module load ucsctools
        bedGraphToBigWig {input.plus} {params.chr_size} {output.plus}
        bedGraphToBigWig {input.minus} {params.chr_size} {output.minus}
        """
