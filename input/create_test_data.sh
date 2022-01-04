module load samtools

ndel1="Chr17:8413131-8490411"
# samtools view -h /home/hsher/293T_rbfox/RBFOX2_293t_GRCh38.1120_2500_CLIP.A04.r1.fq.genome-mappedSo.rmDupSo.merged.r2.bam "chr17:8413131-8490411" | samtools view -Sb - > rbfox2.bam

# samtools view -h /projects/ps-yeolab5/encore/processing/encore_master_hg38/encode4_batch2.4001_CLIP1_CLIP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam "chr17:8413131-8490411" | samtools view -Sb - > nono.bam

#samtools view -h /home/hsher/sara_clip/GSE146878/mapped_reads/CLIP1.Aligned.sortedByCoord.out.bam "chr17:8413131-8490411" | samtools view -Sb - > fmrp.bam

histone="chr6:26055740-26056470"

# turn out to be hg38
slbp_clip_vivo=/home/hsher/f-SHAPE-eCLIP/snake/input/SLBP.CLIP_vivo1.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
samtools index $slbp_clip_vivo
samtools view -h $slbp_clip_vivo  "chr6:26055740-26056470" | samtools view -Sb -> /home/hsher/f-SHAPE-eCLIP/snake/input/slbp_vivo.bam

slbp_clip_vitro=/home/hsher/f-SHAPE-eCLIP/snake/input/SLBP.CLIP_vitro1.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
samtools index $slbp_clip_vitro
samtools view -h $slbp_clip_vitro  "chr6:26055740-26056470" | samtools view -Sb -> /home/hsher/f-SHAPE-eCLIP/snake/input/slbp_vitro.bam

EEF2="chr19:3976056-3985463"
dhx36_n1=/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_N1_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam
dhx36_d1=/home/hsher/rg4_seq/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D1_IP.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam

samtools index $dhx36_n1
samtools view -h $dhx36_n1 "chr19:3976056-3985463" | samtools view -Sb -> /home/hsher/f-SHAPE-eCLIP/snake/input/dhx36_n1.bam

samtools index $dhx36_d1
samtools view -h $dhx36_d1 "chr19:3976056-3985463" | samtools view -Sb -> /home/hsher/f-SHAPE-eCLIP/snake/input/dhx36_d1.bam

