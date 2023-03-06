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
    
    # Write reads to file
    SAM_FILE = pysam.AlignmentFile(sam_path, "r") 
    with open(output_path, 'a') as fa:
        discarded_count = 0
        written_list = []
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
                    discarded_count += 1
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
                # <500 bp alignment in middle of >1kb reads considered a misalignment here
                elif end-start < 500 and len(seq) > 1000:
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
            # Reads above 20000 bp to be filtered for misalignments where "true" alignments are >500 bp
            elif query_length > 20000:
                # Make a list of reads already written, need to do this if >1 alignment of a read is in SAM file
                if name not in written_list:
                    pos_start=[]
                    pos_end=[]
                    # Make a list of "true" alignment positions
                    for pos in algn_dict[name][1:]:
                        if pos[1]-pos[0] > 500:
                            pos_start.append(pos[0])
                            pos_end.append(pos[1])
                    # Read with no "true" alignments
                    if len(pos_start) == 0:
                        fa.write(f">{name}_{1}")
                        fa.write("\n")
                        fa.write(seq)
                        fa.write('\n')
                    if len(pos_start) == 1:
                        if pos_start[0] < 500:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[pos_end[0]:algn_dict[name][0]])
                            fa.write('\n')
                        elif algn_dict[name][0]-pos_end[0] < 500:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[0:pos_start[0]])
                            fa.write('\n')
                        elif pos_end[0]-pos_start[0] < 6500:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[0:pos_start[0]])
                            fa.write('\n')
                            fa.write(f">{name}_{2}")
                            fa.write("\n")
                            fa.write(seq[pos_end[0]:algn_dict[name][0]])
                            fa.write('\n')
                    elif len(pos_start) == 2:
                        if min(pos_start) < 1000:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[min(pos_end):max(pos_start)])
                            fa.write('\n')
                            fa.write(f">{name}_{2}")
                            fa.write("\n")
                            fa.write(seq[max(pos_end):algn_dict[name][0]])
                            fa.write('\n')
                        elif algn_dict[name][0]-max(pos_end) < 1000:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[0:min(pos_start)])
                            fa.write('\n')
                            fa.write(f">{name}_{2}")
                            fa.write("\n")
                            fa.write(seq[min(pos_end):max(pos_start)])
                            fa.write('\n')
                        else:
                            written_list.append(name)
                            fa.write(f">{name}_{1}")
                            fa.write("\n")
                            fa.write(seq[0:min(pos_start)])
                            fa.write('\n')
                            fa.write(f">{name}_{2}")
                            fa.write("\n")
                            fa.write(seq[min(pos_end):max(pos_start)])
                            fa.write('\n')
                            fa.write(f">{name}_{3}")
                            fa.write("\n")
                            fa.write(seq[max(pos_end):algn_dict[i][0]])
                            fa.write('\n') 
            else:
                discarded_count += 1

#print("********************************\n")
#print("\n")
#print(str(discarded_count)+" reads out of "+str(len(algn_dict))+" were discarded")
#print("\n")
#print("********************************\n")

if __name__ == '__main__':
    main()
