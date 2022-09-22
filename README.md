# Hashtag Demultiplex Pipeline

## Files

- **main.df:** Contains Nextflow script, which calls the processes
- **nextflow.config:** Default parameters for the pipeline
- **HTODemux.R:** HTODemux algorithm 
- **HTODemux-visualisation.R:** HTODemux graphs
- **multi-seq.R:** HTODemux algorithm 

### HTODemux 
#### Input

- UMI Matrix.rds File
- Hashtag counts matrix  .rds File. (*)
- Parameters for HTODemux

### MULTI-seq
#### Input
- UMI Matrix.rds File
- Hashtag counts matrix  .rds File. (*)
- Parameters for MULTI-seq

### Hashed Drops
#### Input
- UMI Matrix.rds File
- Parameters for Hashed Drops

### DemuxEM
#### Input
- UMI Matrix hdf5 File
- Hashtag counts matrix - 10x Folder path
- Parameters for DemuxEM

### Hash Solo
#### Input
- Hashtag counts matrix - 10x Folder path
- Parameters for Hash Solo


