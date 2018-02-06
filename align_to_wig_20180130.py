'''
Generates pile-up files (reads counts that map along the genome, with position in nucleotides) wig format. Deals with mismatches (from non-template addition from reverse transcription) and possible multiple alignments, with the explicit purpose of correctly taking into account close paralogs (e.g., EF-Tu). 

Input:
1. An align file (standard bowtie output). Assumes -k at most 2 as a bowtie option.
2. Name of resulting wig files
3. Chromosome name (for wig header)
    
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

# import required packages
import sys, random, copy


def write_wig(inputFile, output_header,chrom_name):
    
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
        
            # if a single alignment, do as usual:
            if (previous_read_name != pprevious_read_name) and (previous_read_name != read_name):
                (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand)
                
                    
            # if two alignments, keep only if a single of the read as a non-template mismatch
            elif previous_read_name == pprevious_read_name:
                                
                # both reads have no mismatches
                if (previous_mismatch_info == '\n') and (pprevious_mismatch_info == '\n'):
                    pass
                # one alignment has a mismatch, but not the other, map the one without mismatch
                elif previous_mismatch_info == '\n':
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand)
                elif pprevious_mismatch_info == '\n':
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,pprevious_mismatch_info,pprevious_fiveprime,pprevious_length,pprevious_strand)
                    
                # non-template addition mismatch
                elif (previous_mismatch_info[0]=='0') and not (pprevious_mismatch_info[0]=='0'):
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,previous_mismatch_info,previous_fiveprime,previous_length,previous_strand)
                elif not (previous_mismatch_info[0]=='0') and (pprevious_mismatch_info[0]=='0'):
                    (forward_5,forward_3,reverse_5,reverse_3) = count_read(forward_5,forward_3,reverse_5,reverse_3,pprevious_mismatch_info,pprevious_fiveprime,pprevious_length,pprevious_strand)
                
                               
    ### Output ###   
    writeOutput_wig(forward_5,output_header+'_5_f.wig',chrom_name)
    writeOutput_wig(forward_3,output_header+'_3_f.wig',chrom_name)
    writeOutput_wig(reverse_5,output_header+'_5_r.wig',chrom_name)
    writeOutput_wig(reverse_3,output_header+'_3_r.wig',chrom_name)
    
    
    
    
def count_read(forward_5,forward_3,reverse_5,reverse_3,mismatch_info,fiveprime,length,strand):
    # count reads (3' and 5' ends of mapped reads) towards the pile-up file.

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
                        
            forward_5[fiveprime+1] = forward_5.get(fiveprime+1,0) + 1
            forward_3[fiveprime+length] = forward_3.get(fiveprime+length,0) + 1        
                    
        elif strand == '-':
            # dealing with non-template addition
            if mismatch == 1:
                length = length - 1
                        
            reverse_5[fiveprime+length]=reverse_5.get(fiveprime+length,0) + 1
            reverse_3[fiveprime+1]=reverse_3.get(fiveprime+1,0) + 1   
        
        
    return (forward_5,forward_3,reverse_5,reverse_3)




def writeOutput_wig(dictionary,File_name,chrom_name):
    # prints the wig file output header. To be changed for species/chromosome name of interest.

    dictionary = dictionary.items()
    dictionary.sort()
    outFile = open(File_name, 'w')
    outFile.write('track type=wiggle_0\n')
    outFile.write('variableStep chrom='+chrom_name+'\n')
    
    #if species=='bsub':
    #    outFile.write('variableStep chrom=NC_000964.3\n')
    #elif species=='ecoli':
    #    outFile.write('variableStep chrom=NC_000913_2\n')
    #elif species=='vnat1':
    #    outFile.write('variableStep chrom=CP009977.1\n')
    #elif species=='vnat2':
    #    outFile.write('variableStep chrom=CP009978.1\n')
    #elif species=='caulobacter':
    #    outFile.write('variableStep chrom=CP001340.1\n')
        
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
    
    # chromosome name for wig file header
    chrom_name = sys.argv[3]
       
    ### generating the .wig files
    write_wig(alignFile, output_header,chrom_name)
    
    

    

