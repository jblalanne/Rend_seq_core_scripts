'''
Generates pile-up files (reads counts that map along the genome, with position in nucleotides) wig format, correcting for peak shadows from identified peaks in the raw data. Deals with mismatches (from non-template addition from reverse transcription) and possible multiple alignments, with the explicit purpose of correctly taking into account close paralogs (e.g., EF-Tu). 

Input:
1. An align file (standard bowtie output). Assumes -k at most 2 as a bowtie option.
2. Name of resulting wig files.
3. end-enrichement file (three columns: 1 position, 2 end-enrichment, 3 type of ends (1:3'f, 2:3'r, 3:5'f, 4:5'r).
4. Chromosome name (for wig header).
    
Output: wig files

The structure of the source file (bowtie default) contains 11 columns of information.
    0. The read name, which is completely unique.
    1. The RefSeq strand of alignment, which is either '+' for forward or '-' for reverse.
    2. The RefSeq chromosome name.
    3. The RefSeq position of the 5' end.
    4. The sequencing read in nucleotides.
    5. The sequencing quality score.
    6. number of other instances the same sequence aligned against the same reference characters (generally not the number of other alignments). 
    7. Mismatches in the alignment.

'''

# importing required packages
import sys, random, copy
import numpy as np


def write_wig_no_shadow(inputFile, output_header,enrichment_data,chrom_name):
    # parses the bowtie align file and count reads.
    
    
    ### Set initial variables ###   
    forward_5, forward_3 = {}, {}
    reverse_5, reverse_3 = {}, {}
    
    ### Run file ###
    inFile = open(inputFile, 'r')
    n_line = 0

    for line in inFile:
        n_line = n_line + 1
        
        # storing previous lines to assess number of alignment
        if n_line>2:
    	    pprevious_length = previous_length
    	    pprevious_fiveprime = previous_fiveprime
    	    pprevious_strand = previous_strand
    	    pprevious_n_other_alignments = previous_n_other_alignments
    	    pprevious_read_name = previous_read_name
    	    pprevious_mismatch_info = previous_mismatch_info
    	    
        
        if n_line>1:
    	    previous_length = length
    	    previous_fiveprime = fiveprime
    	    previous_strand = strand
    	    previous_n_other_alignments = n_other_alignments
    	    previous_read_name = read_name
    	    previous_mismatch_info = mismatch_info
    	    
    	
    	# parsing each line of the bowtie output.
        fields = line.split('\t')
        read_name = fields[0]
        length = len(fields[4])
        fiveprime = int(fields[3])
        strand = str(fields[1])
        n_other_alignments = fields[6]
        mismatch_info = fields[7]
        
        
        if n_line>2:
        
            # if a single alignment, count read directly:
            if (previous_read_name != pprevious_read_name) and (previous_read_name != read_name):
                (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand,enrichment_data)
                
                    
            # if two alignments, keep only if a single of the read as a non-template mismatch
            elif previous_read_name == pprevious_read_name:
                                
                # both reads have no mismatches
                if (previous_mismatch_info == '\n') and (pprevious_mismatch_info == '\n'):
                    pass
                # one alignment has a mismatch, but not the other, map the one without mismatch
                elif previous_mismatch_info == '\n':
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand,enrichment_data)
                elif pprevious_mismatch_info == '\n':
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,pprevious_mismatch_info,pprevious_fiveprime,pprevious_length,pprevious_strand,enrichment_data)
                    
                # non-template addition mismatch
                elif (previous_mismatch_info[0]=='0') and not (pprevious_mismatch_info[0]=='0'):
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand,enrichment_data)
                elif not (previous_mismatch_info[0]=='0') and (pprevious_mismatch_info[0]=='0'):
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,pprevious_mismatch_info,pprevious_fiveprime,pprevious_length,pprevious_strand,enrichment_data)
                
                               
    ### Done reading align file, printing output ###   
    writeOutput_wig(forward_5,output_header+'_5_f.wig',chrom_name)
    writeOutput_wig(forward_3,output_header+'_3_f.wig',chrom_name)
    writeOutput_wig(reverse_5,output_header+'_5_r.wig',chrom_name)
    writeOutput_wig(reverse_3,output_header+'_3_r.wig',chrom_name)
    
    
    
    
def count_read(forward_5,forward_3,reverse_5,reverse_3,mismatch_info,fiveprime,length,strand,enrichment_data,distance_buffer=3):
    # count reads (3' and 5' ends of mapped reads) towards the pile-up file, down-weighting positions identified as peaks.

    # determining the nature of mismatch. Only keepin mismatches resulting (likely) from non-template addition of reverse transcription (mismatch at 5' end of cDNA) or bad read position (N). 
    if mismatch_info == '\n':
        mismatch = 0
    elif '0' == mismatch_info[0]:
        mismatch = 1
    elif '>N' in mismatch_info:
        mismatch = 0
    else:
        mismatch = 2
    
    # only keeping reads in desired range and with mismatch type (above).  
    if (14 < length < 45) and (mismatch < 2):
    
        if strand == '+':
        
            # dealing with non-template addition
            if mismatch == 1:
                fiveprime += 1
                length = length - 1
                
            position_5 = fiveprime+1
            position_3 = fiveprime+length
            
            # is the position a peak (within distance_buffer-1 nt)? Look through end-enrichment file.
            bool_5 = np.where( (abs(enrichment_data[:,1]-position_5)<distance_buffer) & (enrichment_data[:,0]==3))
            bool_3 = np.where( (abs(enrichment_data[:,1]-position_3)<distance_buffer) & (enrichment_data[:,0]==1))
    
            if len(bool_5[0])==0 and len(bool_3[0])==0:			# not a peak, add 1 to both 5' and 3' ends
                forward_5[position_5] = forward_5.get(position_5,0) + 1.0
                forward_3[position_3] = forward_3.get(position_3,0) + 1.0
            elif len(bool_5[0])==1 and len(bool_3[0])==0:		# 5' peak, add 1 to 5' end, 1/end-enrichment at 3' end
                forward_5[position_5] = forward_5.get(position_5,0) + 1.0
                forward_3[position_3] = forward_3.get(position_3,0) + 1.0/enrichment_data[bool_5[0][0],2]
            elif len(bool_5[0])==0 and len(bool_3[0])==1:		# 3' peak, add 1 to 3' end, 1/end-enrichment at 5' end
                forward_5[position_5] = forward_5.get(position_5,0) + 1.0/enrichment_data[bool_3[0][0],2]
                forward_3[position_3] = forward_3.get(position_3,0) + 1.0
            elif len(bool_5[0])==1 and len(bool_3[0])==1:		# double peak, likely PCR jackpot, 1/end-enrichment for both 5' and 3' ends
                forward_5[position_5] = forward_5.get(position_5,0) + 1.0/enrichment_data[bool_3[0][0],2]
                forward_3[position_3] = forward_3.get(position_3,0) + 1.0/enrichment_data[bool_5[0][0],2]
                      
                    
        elif strand == '-':
        
            # dealing with non-template addition
            if mismatch == 1:
                length = length - 1
                
            position_3 = fiveprime+1
            position_5 = fiveprime+length  
            
            # is the position a peak (within distance_buffer-1 nt)? Look through end-enrichment file.
            bool_5 = np.where( (abs(enrichment_data[:,1]-position_5)<distance_buffer) & (enrichment_data[:,0]==4))
            bool_3 = np.where( (abs(enrichment_data[:,1]-position_3)<distance_buffer) & (enrichment_data[:,0]==2))
    
            if len(bool_5[0])==0 and len(bool_3[0])==0:			# not a peak, add 1 to both 5' and 3' ends
                reverse_5[position_5] = reverse_5.get(position_5,0) + 1.0
                reverse_3[position_3] = reverse_3.get(position_3,0) + 1.0
            elif len(bool_5[0])==1 and len(bool_3[0])==0:		# 5' peak, add 1 to 5' end, 1/end-enrichment at 3' end
                reverse_5[position_5] = reverse_5.get(position_5,0) + 1.0
                reverse_3[position_3] = reverse_3.get(position_3,0) + 1.0/enrichment_data[bool_5[0][0],2]
            elif len(bool_5[0])==0 and len(bool_3[0])==1:		# 3' peak, add 1 to 3' end, 1/end-enrichment at 5' end
                reverse_5[position_5] = reverse_5.get(position_5,0) + 1.0/enrichment_data[bool_3[0][0],2]
                reverse_3[position_3] = reverse_3.get(position_3,0) + 1.0
            elif len(bool_5[0])==1 and len(bool_3[0])==1:		# double peak, likely PCR jackpot, 1/end-enrichment for both 5' and 3' ends
                reverse_5[position_5] = reverse_5.get(position_5,0) + 1.0/enrichment_data[bool_3[0][0],2]
                reverse_3[position_3] = reverse_3.get(position_3,0) + 1.0/enrichment_data[bool_5[0][0],2]
        
    return (forward_5,forward_3,reverse_5,reverse_3)


def writeOutput_wig(dictionary,File_name,chrom_name):
    # prints the wig file output header. To be changed for species/chromosome name of interest.

    dictionary = dictionary.items()
    dictionary.sort()
    outFile = open(File_name, 'w')
    outFile.write('track type=wiggle_0\n')
    outFile.write('variableStep chrom='+chrom_name+'\n')
        
    for y in dictionary:
        position = str(y[0])
        read_counts = str(y[1])
        outFile.write(position + '\t' + read_counts + '\n')
    outFile.close()

    
        
if __name__ == '__main__':
    
    ### input from command call
    
    # name of .align file
    alignFile = sys.argv[1]
    
    # name of resulting .wig files
    output_header = sys.argv[2]
    
    # name of enrichment factor file
    enrichment_file = sys.argv[3]
    
    # chromosome name for wig file header
    chrom_name = sys.argv[4]
    
    ### starting script
    
    # read in the enrichment factor file
    fid = open(enrichment_file) 
    enrichment_data = np.loadtxt(enrichment_file, delimiter="\t")
    
    # generating the .wig files
    write_wig_no_shadow(alignFile, output_header,enrichment_data,chrom_name)
    
    
    
    

