
# This pipeline is now deprecated, please use [NanoSim](https://github.com/bcgsc/NanoSim) to simulate both genomic and transcriptomic ONT reads!

# Trans-NanoSim
Oxford nanopore transcriptome read simulator

Trans-NanoSim is a fast and scalable read simulator that captures the technology-specific features of ONT transcriptome reads (cDNA / directRNA), and allows for adjustments upon improvement of nanopore sequencing technology.  

## Dependencies
minimap2 (Tested with version 2.10)  
LAST (Tested with version 581 and 916)  
GenomeTools : http://genometools.org/  
Python (2.7 or >= 3.4)  
Python packages:  
* six  
* numpy (Tested with version 1.10.1 or above)
* HTSeq  
* scipy (Tested with verson 1.0.0)

## Usage
Trans-NanoSim is implemented using Python for error model fitting, read length analysis, and simulation. The first step of Trans-NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference transcriptome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

### 1. Characterization stage  
Characterization stage takes a reference transcriptome / reference genome and a training read set in FASTA format as input and aligns these reads to the references using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats.  

__Usage:__  
```
./read_analysis.py <options>  
    [options]:  
    -h : print usage message  
    -i : training ONT reads (cDNA / directRNA), must be fasta files  
    -rg : reference genome of the training reads  
    -rt : reference transcriptome for the training reads
    -annot : annotation file in ensembl GTF/GFF3 formats (required for Intron Retention modeling)  
    -a : Aligner to be used: minimap2 or LAST, default = 'minimap2' 
    -ga : genome alignment file in sam or maf format (optional)
    -ta : transcriptome alignment file in sam or maf format (optional)  
    -o : The prefix of output files, default = 'training'  
    --no_model_fit : disable the modeling fitting step  
    --no_intron_retention : disable the Intron Retention modeling step  
    --detect_IR : detect Intron Retention events and outputs the markov model for it  
    --quantify : quantify transcripts expression profile
    -
```

### 2. Simulation stage  
Simulation stage takes reference genome and read profiles as input and outputs simulated reads in FASTA fomat.  

__Usage:__  
```
./simulator.py [command] <options>  
   [command]:  
    circular | linear  
    # Do not choose 'circular' when there is more than one sequence in the reference  
    <options>:  
    -h : print usage message
    -rt : reference transcriptome
    -rg : reference genome
    -e : expression profile in the specified format  
    -c : the prefix and address of training set profiles, same as the output prefix in read_analysis.py, default = training  
    -o : The prefix and address of output files, default = 'simulated'  
    -n : Number of generated reads, default = 20,000 reads  
    --perfect: Output perfect reads, no mutations, default = False  
    --KmerBias: prohibits homopolymers with length >= 6 bases in output reads, can be omitted  
```


## Explaination of output files  
### 1. Characterization stage
1. `training_aligned_length_ecdf` Length distribution of aligned regions on aligned reads  
2. `training_aligned_reads_ecdf` Length distribution of aligned reads  
3. `training_align_ratio` Empirical distribution of align ratio of each read  
4. `training_besthit.maf` The best alignment of each read based on length  
5. `training_match.hist/training_mis.hist/training_del.hist/training_ins.hist` Histogram of match, mismatch, and indels  
6. `training_first_match.hist` Histogram of the first match length of each alignment  
7. `training_error_markov_model` Markov model of error types  
8. `training_ht_ratio` Empirical distribution of the head region vs total unaligned region  
9. `training.maf` The output of LAST, alignment file in MAF format  
10. `training_match_markov_model` Markov model of the length of matches (stretches of correct base calls)  
11. `training_model_profile` Fitted model for errors  
12. `training_processed.maf` A re-formatted MAF file for user-provided alignment file  
13. `training_unaligned_length_ecdf` Length distribution of unaligned reads  
14. `training_error_rate.tsv` Mismatch rate, insertion rate and deletion rate

### 2. Simulation stage  
1. `simulated.log`  
  Log file for simulation process  
  
2. `simulated_reads.fasta`  
  FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
  
3. `simulated_error_profile`  
  Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases. 
