import pandas as pd
try:
    MANIFEST=config['MANIFEST']
    manifest = pd.read_table(MANIFEST, index_col = False, sep = ',')
    print(manifest.head())
    sample_labels = manifest['Sample'].tolist()
    print(sample_labels)
    
    workdir: config['WORKDIR']

    print(os.getcwd())

    os.mkdir('error_files')
except Exception as e:
    print(e)

# snakemake -s SnakeMain.py --configfile config/encode.yaml -j 30 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
# snakemake -s SnakeMain.py --configfile config/encode.yaml -j 30 --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"
module tracks:
    snakefile:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        '/home/hsher/projects/make_tracks/Snake_CLIP_track.py'
    config:
        config

rule all:
    input:
        expand("CITS/{sample_label}.pos.bw", sample_label = sample_labels)+
        #expand("CIMS/{sample_label}.pos.bw", sample_label = sample_labels)+
        expand("coverage/{sample_label}.pos.bw", sample_label = sample_labels)
    
use rule* from tracks as track_*
