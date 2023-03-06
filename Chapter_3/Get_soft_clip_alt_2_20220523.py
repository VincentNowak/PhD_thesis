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



    # clean up if out file alread exists 
    if output_path.is_file:
        output_path.unlink(missing_ok=True)

    # Make dictionary with alignment information per read as reference
    # read qualities (as previously) but also CIGAR strings are disregarded in this version 
    SAM_FILE = pysam.AlignmentFile(sam_path, "r")
    algn_dict = {}
    # Note that this loops through alignments and not reads!
    for read in SAM_FILE:
        if read.query_name in algn_dict:
            algn_dict[read.query_name].append((read.query_alignment_start,read.query_alignment_end))
        else:
            algn_dict[read.query_name] = [read.infer_query_length(),(read.query_alignment_start,read.query_alignment_end)]
    
    SAM_FILE = pysam.AlignmentFile(sam_path, "r") 
    # Write reads to file
    with open(output_path, 'a') as fa:
        for read in SAM_FILE:
            seq = read.query_sequence
            name = read.query_name
            query_length = algn_dict[read.query_name][0]
            # Unmapped reads
            if query_length == None:
                fa.write(f">{name}_{1}")
                fa.write("\n")
                fa.write(seq)
                fa.write('\n')
           # Single alignment
            elif len(algn_dict[name]) == 2:
                start = algn_dict[name][1][0]
                end = algn_dict[name][1][1]
                if len(seq) == end-start:
                    continue
                elif start < 500 and len(seq[end:query_length]) > 1000:
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq[algn_dict[name][1][1]:algn_dict[name][0]])
                    fa.write('\n')
                elif query_length-end < 500 and len(seq[0:start]) > 1000:
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq[0:algn_dict[name][1][0]])
                    fa.write('\n')
                # <250 bp alignment in middle of >1kb reads considered a misalignment here
                elif end-start < 250 and len(seq) > 1000:
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq)
                    fa.write('\n')
                # Using the 3500 from below doesn't change number much
                elif end-start < 6500 and len(seq[0:start]) > 1000 and len(seq[end:query_length]) > 1000:                
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq[0:algn_dict[name][1][0]])
                    fa.write('\n')
                    fa.write(f">{name}_{2}")
                    fa.write("\n")
                    fa.write(seq[algn_dict[name][1][1]:algn_dict[name][0]])
                    fa.write('\n')
            # Two alignments and >20 kb (mostly insert; changing to >10kb doesn't change numbers much)
            elif len(algn_dict[name]) == 3 and query_length > 20000:
                # determine max and min values in alignment
                max_left = max(algn_dict[name][1:],key=lambda item:item[1])[0]
                min_left = min(algn_dict[name][1:],key=lambda item:item[1])[0]
                min_right = min(algn_dict[name][1:],key=lambda item:item[1])[1]
                max_right = max(algn_dict[name][1:],key=lambda item:item[1])[1]
                # case where vector is at both ends and full insert is sequenced
                if max_left-min_right > 20000:
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq[min_right:max_left])
                    fa.write('\n')
                # case where vector is in middle
                if max_right-min_left < 6500:
                    fa.write(f">{name}_{1}")
                    fa.write("\n")
                    fa.write(seq[0:min_left])
                    fa.write('\n')
                    fa.write(f">{name}_{2}")
                    fa.write("\n")
                    fa.write(seq[max_right:query_length])
                    fa.write('\n')


if __name__ == '__main__':
    main()
