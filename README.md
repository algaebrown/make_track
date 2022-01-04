# make_track
Making eCLIP coverage, CITS, CIMS tracks

example data are in the `input/` folder.

`snakemake -j 30 -s Snake_CLIP_track.py --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn=1 -q home-yeo -e {params.error_out_file} -o /dev/null"` to run the pipeline
