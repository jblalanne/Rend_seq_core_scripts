# Rend-seq core scripts

Core scripts used for analysis of Rend-seq data.
Included MATLAB and python scripts allow for the generation of pile-up (.wig) files from aligned files (bowtie standard output), the removal of peak shadows, the identification of peaks, and the display and quantification of mRNA isoforms from a region of interest. 

## wig file generation
Detailed below are the steps to generate wig files for Rend-seq data, starting from a standard bowtie align file.

### Raw wig file
To obtain the raw Rend-seq wig, run:
```
python align_to_wig_20180130.py my_align_file.align my_wigs my_wig_chrom_name
```
The above command will generate four wig files: my_wigs_3f.wig, my_wigs_3r.wig, my_wigs_5f.wig, and my_wigs_5r.wig. These respectively correspond to the read counts for the 3' end of mapped reads on the forward strand, 3' end of mapped reads on the reverse strand, 5' end of mapped reads on the forward strand, and 5' end of mapped reads on the reverse strand. Each wig file will have as header
```
track type=wiggle_0
variableStep chrom=my_wig_chrom_name
```
In the above command, my_align_file.align is the standard bowtie align file (typically with -v 1 -k 1, or -v 1 -m 2 -k 2 options, the above script assumes at most two alignments are reported, or option -k is less than or equal to 2). Refer to the [bowtie manual](http://bowtie-bio.sourceforge.net/manual.shtml) for details.  

The generated wig files can be used for data visualization in a genome browser, such as [MochiView](http://www.johnsonlab.ucsf.edu/mochi/).

### End-enrichment file for peak shadow removal
The raw wig files generated above can be used to identify peaks in the reads counts, corresponding to ends of transcripts. The signal enrichment at these ends can then be computed and used to remove peak shadows (see Methods of corresponding publication for details). 

The MATLAB script `Rend_seq_end_enrichment_file_header_20180206.m` contains an example (using Rend-seq data from exponentially growing *Bacillus subtilis* in LB, whose raw wigs are included in this repository). To use the MATLAB scripts, download all MATLAB files from this submission and add the directory (and sub-directories) containing these files to the MATLAB path (or run the script from the directory containing these files). You will have to change variable `data_dir` to the full path to the directory containing the wig files. 

`Rend_seq_end_enrichment_file_header_20180206.m` first imports the data in MATLAB, then computes the peak z score and identify peaks, and then generates the end-enrichment file (which contains peak type (3f, 3r, 5f or 5r), position, and end-enrichment). For example included, the file generated will be `Bacillus_subtilis_WT_Rend_seq_LB_25s_frag_pooled_end_enrichment_shadow_removal.txt`. This file can then be used to generate the shadow removed wig. See script comments and Methods of the associated publication for details.

Any peak missed by the automated identification (e.g., if a 3' end has 3 to 5 different possible neighboring positions) will not have its shadow removed. To also remove the shadow of a spuriously missed peak (by the automated analysis), information can be added to the end-enrichment file manually. 

### Shadow removed wig file

The end-enrichment file (see above) and the original align file can then be used to generate shadow removed wig files. To do so, run:
```
python align_to_wig_no_shadow_20180130.py my_align_file.align my_wigs_no_shadow my_file_end_enrichment_shadow_removal.txt my_wig_chrom_name
```
This is similar to the script used for the raw wig file generation, except that now the end-enrichment file (described above) is also required. This script generates four wig files (similar as before): my_wigs_no_shadow_3f.wig, my_wigs_no_shadow_3r.wig, my_wigs_no_shadow_5f.wig, and my_wigs_no_shadow_5r.wig. 


## mRNA isoforms display and quantification

mRNA isoforms from complex operons can be identified and quantified from the shadow removed wigs (the scripts also works with raw wig files). The MATLAB script `Rend_seq_isoform_header_20180130.m` contains working examples and can easily be modified to display any region of interest. It uses shadow removed wig files from Bacillus subtilis included in the current repository. It also requires a GenBank annotation file (obtained for a genome of interest by downloading coding sequences in the FASTA Protein format). That of [*Bacillus subtilis*](https://www.ncbi.nlm.nih.gov/nuccore/NC_000964) (`NC_000964.3.faa`) is included in the repository. 

As before, the scripts must be downloaded and the containing directory (and sub-directories) must be added to the MATLAB path.  The variable `data_dir` and `annotation_dir` must then be changed to the directories containing the shadow removed wig files and GenBank annotation files respectively. 

The script first identifies peaks (corresponding to ends of transcripts) in the region of interest, and then estimates the isoform relative abundances as described in the Methods of the associated publication. The current example is for the *rpsP* operon. The analysis can be performed for other regions by changing variables `start_region`, `stop_region`, and `strand_region`. See further description, examples, and comments in the script. 

