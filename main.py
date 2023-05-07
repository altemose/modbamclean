import numpy as np
import pysam
import argparse
from array import array

# get filepaths from command line
parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="bam file", default="test.bam", required=False)
parser.add_argument("--mapq", help="mapping quality threshold (int 0-60)", default=10, required=False)
parser.add_argument("--baseq", help="basecall quality threshold (int 0-60)", default=8, required=False)
parser.add_argument("--modq", help="modification score threshold (int 0-255)", default=155, required=False)
parser.add_argument("--out", help="output bam file", default="output.bam", required=False)
args = parser.parse_args()

bam = args.bam
output = args.out
mapq_score = args.mapq
baseq_score = args.baseq
mod_score = args.modq
np.set_printoptions(threshold=np.inf)

input_bam_file = pysam.AlignmentFile(
    bam, "rb")
output_bam_file = pysam.AlignmentFile(
    output,
    "wb", template=input_bam_file, header=input_bam_file.header)

# iterate through each alignment in the BAM file
num = 0
for read in input_bam_file.fetch():
    # mapq filter
    if read.is_unmapped or read.is_secondary or read.is_supplementary or (read.mapping_quality < mapq_score):
        continue

    else:
        if(read.has_tag('Ml')==False):
            output_bam_file.write(read)
            continue

        # baseq filter
        read_sequence=read.get_forward_sequence()
        baseq_pass = np.array(read.get_forward_qualities()).__ge__(baseq_score) #get boolean array for positions with baseq >= thresh
        sequence_length = baseq_pass.size
        sequence_length2 = len(read_sequence)

        # ref=A filter
        aligned_pairs = np.array(read.get_aligned_pairs(with_seq=True))
        aligned_pairs = aligned_pairs[aligned_pairs[:, 0].__ne__(None)] # eliminate gaps in the aligned read sequence
        refA_pass = aligned_pairs[:, 2].__eq__('A') #get boolean array for positions with A in ref

        # mod score filter
        mod = np.array(read.modified_bases_forward[('A', 0, 'a')]) #gets A mod scores
        mod_thresh_pass = mod[mod[:,1].__ge__(mod_score),:] #ignores those <thresh
        if(mod_thresh_pass.size==0):
            read.set_tag('Mm',value=None)
            read.set_tag('Ml',value=None)
            output_bam_file.write(read)
            continue
        mod_scores = np.zeros(sequence_length) #initialize array of zeros
        mod_scores[mod_thresh_pass[:,0]] = mod_thresh_pass[:,1] #fill in positions above thresh only

        #apply baseq and refA filters to mod_scores already passing modq threshold
        mod_scores[np.logical_and(refA_pass, baseq_pass)==False] = 0
        if (sum(mod_scores) == 0):
            read.set_tag('Mm',value=None)
            read.set_tag('Ml',value=None)
            output_bam_file.write(read)
            continue

        #create array of final passing positions and scores, reformat to mm and ml
        mod_positions=np.argwhere(mod_scores > 0)
        mltag=array('B',mod_scores[mod_positions].astype(int)[:,0])#coerce to type array.array for pysam set_tag compatibility

        # Find indices of all A's in sequence with numpy (Reinterpret str as a char buffer for speed)
        A_positions = np.argwhere(
            np.frombuffer(
                read_sequence.encode(), dtype=np.uint8) == ord('A'))

        # create the mmtag using numpy searchsorted and ediff1d functions
        mmtag='A+a,'+','.join((np.ediff1d(
            np.insert(
                np.searchsorted(
                    A_positions[:, 0], mod_positions[:, 0], side='left'),0,0)) - 1).astype(str))
        #print(mmtag)

        read.set_tag('Mm',value=mmtag)  # store in read object
        read.set_tag('Ml',value=mltag)  # store in read object

        output_bam_file.write(read)

    #num += 1
    #if num > 1:
    #    break


input_bam_file.close()
output_bam_file.close()
