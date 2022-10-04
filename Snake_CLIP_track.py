#configfile: "config.yaml"
# snakemake -j 30 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
# snakemake -j 200 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn=1 -q home-yeo" --configfile config/downsample_bam.yaml --keep-going 
# snakemake -j 20 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time}:00:00 -l nodes=1:ppn=1 -q home-yeo" --configfile config/full_bam.yaml --keep-going 


import pandas as pd
import glob




# this might not work for single-end?
# yes when single end we only have header
rule extract_read_two:
    input:
        #bam="input/{sample_label}.bam"
        bam=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["bam"].values[0]),
    output:
        read2="processed_bam/{sample_label}.r2.bam",
        read1="processed_bam/{sample_label}.r1.bam"
    params:
        run_time="2:00:00",
        cores = 1,
        error_out_file = "error_files/extract_read2",
        out_file = "stdout/extract_read2.{sample_label}"
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
        pos_bam="processed_bam/{sample_label}.r2.pos.bam",
        neg_bam="processed_bam/{sample_label}.r2.neg.bam"
    params:
        run_time="2:00:00",
        cores = 1,
        error_out_file = "error_files/strand.bam.{sample_label}",
        out_file = "stdout/strand.bam.{sample_label}"
    shell:
        """
        # make strand specific
        module load samtools
        
        samtools view -h -F 0x10 {input.bam} | samtools view -Sb - > {output.pos_bam}
        samtools view -h -f 0x10 {input.bam}| samtools view -Sb - > {output.neg_bam}
        """


rule bam_pileup:
    input:
        pos_bam="processed_bam/{sample_label}.r2.pos.bam",
        neg_bam="processed_bam/{sample_label}.r2.neg.bam"
    output:
        pos_pileup="CIMS/{sample_label}.pos.pileup",
        neg_pileup="CIMS/{sample_label}.neg.pileup",
        
    params:
        run_time="2:00:00",
        genomefa=config['GENOME_FA'],
        error_out_file = "error_files/bam_pileup",
        out_file = "stdout/pileup.{sample_label}",
        cores = 1,
    shell:
        """
        module load samtools
        
        samtools mpileup -s  {input.pos_bam} -f {params.genomefa} > {output.pos_pileup}
        samtools mpileup -s  {input.neg_bam} -f {params.genomefa} > {output.neg_pileup}
        module unload samtools
        """
rule CIMS_bw:
    input:
        pos_pileup="CIMS/{sample_label}.pos.pileup",
        neg_pileup="CIMS/{sample_label}.neg.pileup",
    output:
        pos="CIMS/{sample_label}.pos.bw",
        neg="CIMS/{sample_label}.neg.bw"
        
    conda:
        "/envs/metadensity.yaml"
    params:
        run_time="2:00:00",
        cores = 1,
        chr_size=config['CHROM_SIZES'],
        error_out_file = "error_files/CIMS_bw",
        out_file = "stdout/CIMS_bw.{sample_label}"
    shell:
        """
        # coverage
        python pileupToMismatchBw.py {input.pos_pileup} {params.chr_size} {output.pos}
        python pileupToMismatchBw.py {input.neg_pileup} {params.chr_size} {output.neg}
        """
rule CITS_bam_to_bedgrah:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        pos="CITS/{sample_label}.pos.bedgraph",
        neg="CITS/{sample_label}.neg.bedgraph"
    params:
        run_time="2:00:00",
        error_out_file = "error_files/CITS_bedgraph",
        cores = 1,
        out_file = "stdout/CITS_bedgraph.{sample_label}"
    shell:
        """
        module load bedtools;
        set +o pipefail;
        
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + -5 > {output.pos}
        bedSort {output.pos} {output.pos}

        bedtools genomecov -ibam {input.bam} -bg -scale -1 -strand - -5 > {output.neg}
        bedSort {output.neg} {output.neg}
        """
rule CITS_bedgraph_to_bw:
    input:
        pos="CITS/{sample_label}.pos.bedgraph",
        neg="CITS/{sample_label}.neg.bedgraph"
    output:
        pos="CITS/{sample_label}.pos.bw",
        neg="CITS/{sample_label}.neg.bw"
    params:
        run_time="2:00:00",
        chr_size=config['CHROM_SIZES'],
        error_out_file = "error_files/CITS_bw",
        out_file = "stdout/CITS_bw.{sample_label}",
        cores = 1,
    shell:
        """
        module load ucsctools
        # filter ERCC
        cat {input.pos} | grep -v ERCC > {input.pos}.filtered
        cat {input.neg} | grep -v ERCC > {input.neg}.filtered

        bedGraphToBigWig {input.pos}.filtered {params.chr_size} {output.pos}
        bedGraphToBigWig {input.neg}.filtered {params.chr_size} {output.neg}
        """
rule COV_bam_to_bedgraph:
    input:
        bam="processed_bam/{sample_label}.r2.bam"
    output:
        pos="coverage/{sample_label}.pos.bedgraph",
        neg="coverage/{sample_label}.neg.bedgraph"
    params:
        run_time="1:00:00",
        error_out_file = "error_files/coverage_bedgraph",
        out_file = "stdout/COV_bedgraph.{sample_label}",
        cores = 1,
    shell:
        """
        module load bedtools;
        set +o pipefail;

        # coverage
        bedtools genomecov -ibam {input.bam} -bg -scale 1 -strand + -split > {output.pos}
        bedSort {output.pos} {output.pos}
        bedtools genomecov -ibam {input.bam} -bg -scale -1 -strand - -split > {output.neg}
        bedSort {output.neg} {output.neg}
        """
rule COV_bedgraph_to_bw:
    input:
        pos="coverage/{sample_label}.pos.bedgraph",
        neg="coverage/{sample_label}.neg.bedgraph"
    output:
        pos="coverage/{sample_label}.pos.bw",
        neg="coverage/{sample_label}.neg.bw"
    params:
        run_time="2:00:00",
        chr_size=config['CHROM_SIZES'],
        error_out_file = "error_files/CITS_bw",
        out_file = "stdout/COV_bw.{sample_label}",
        cores = 1,
    shell:
        """
        module load ucsctools
        # filter ERCC
        cat {input.pos} | grep -v ERCC > {input.pos}.filtered
        cat {input.neg} | grep -v ERCC > {input.neg}.filtered

        bedGraphToBigWig {input.pos}.filtered {params.chr_size} {output.pos}
        bedGraphToBigWig {input.neg}.filtered {params.chr_size} {output.neg}
        """
