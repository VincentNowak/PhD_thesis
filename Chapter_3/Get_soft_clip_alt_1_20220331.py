from pathlib import Path
import pysam
import argparse

'''feed in a SAM file from minimap2 (reads mapped to vector).
Will write out the softclipped regions, which should be the insert seqs.
Positional args will be the path the .sam file and the output prefix
'''

def main():
    
    arg_parsers = argparse.ArgumentParser(description='cut out vector seqs and return inserts')

    arg_parsers.add_argument('-s',
                           type=str,
                           help='the path to sam file')

    arg_parsers.add_argument('-o',
                           type=str,
                           help='output')

    args = arg_parsers.parse_args()

    sam_path = Path(args.s)

    output_path = Path(args.o)




    def get_regions(query_sequence:str, cigar_string:str) -> dict:
        '''move across the read and bin by alignment type (clipped or aligned)''' 

        read_consuming_op_ints = [0, 1, 4, 7, 8]
        read_consuming_ops = ["M", "I", "S", "=", "X"]
        # to hold clipped and aligned regions
        seq_extract = {op: [] for op in read_consuming_ops}
        # current position in read walk
        pos = 0 
        for c in cigar_string:
            op_int, length = c

            if op_int in read_consuming_op_ints:

                # get the op as a str char 
                op = read_consuming_ops[read_consuming_op_ints.index(op_int)]
                if query_sequence:
                    seq_region = query_sequence[pos:pos+length]
                    seq_extract[op].append(seq_region)
                pos += length

        return seq_extract
    
    
    # open sam file
    SAM_FILE = pysam.AlignmentFile(sam_path, "r") # open samfile

    # clean up if out file alread exists 
    if output_path.is_file:
        output_path.unlink(missing_ok=True)

    # open and write to output file
    with open(output_path, 'a') as fa:
        
        #count the reads for debug
        read_number = 0
        for read in SAM_FILE:# open samfile
            # get the mapping data for each read 
            cigar_string = read.cigar
            query_sequence = read.query_sequence
            query_name = read.query_name

            # This was added to also write the unmapped reads (without any vector) to the fasta file
            if read.is_unmapped == True:
                fa.write(f">{query_name}")
                fa.write("\n")
                fa.write(query_sequence)
                fa.write('\n')

            if query_sequence and cigar_string:
                seq_extract_dict = get_regions(query_sequence, cigar_string)
                #print(read_number)
                #print(query_name)
                #print(cigar_string)
                #print("")
                read_number += 1

                # grab the softclipped bits and write them to the out file 
                for i, soft_clip in enumerate(seq_extract_dict["S"]):
                    fa.write(f">{query_name}_{i}")
                    fa.write("\n")
                    fa.write(soft_clip)
                    fa.write('\n')
                
        print("finished")
        
if __name__ == '__main__':
    main()
