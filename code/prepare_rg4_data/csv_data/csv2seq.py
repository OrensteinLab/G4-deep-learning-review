from extract_most_prominent_transcript import get_transcript_and_position, transcript
import sys
import os
import pandas as pd
import pickle
from Bio import SeqIO
import re

'''
@author: Maor Turner
'''

# Execution:
# under csv_data directory run:
# csv2seq.py  <chr_dict_path> <transcripts_fasta_path>

FLANK_SIZE = 47
SEQ_SIZE = 30

# get path from command line
if len(sys.argv) != 3:
    print("ERROR:\n"
          "Execution: csv2seq.py  <chr_dict_path> <transcripts_fasta_path>\n"
          "Run under csv_data directory")
    exit(1)
chr_dict_path = sys.argv[1]
transcripts_fasta_path = sys.argv[2]

# set output directory
output_dir = "/out/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# load transcripts locations and tsl data
with open(chr_dict_path, 'rb') as f:
    chr_dict = pickle.load(f)

# load transcripts in to dictionary and set keys
original_transcripts_dict = SeqIO.to_dict(SeqIO.parse(transcripts_fasta_path, "fasta"))
transcripts_dict = {}
for key in original_transcripts_dict:
    m = re.search(r'(ENST\d+\.\d*)|', key)
    assert m, f"Did not find transcript ID in {key}"
    transcripts_dict[m.group(1)] = original_transcripts_dict[key]
del original_transcripts_dict

# you can also loop over all 3 data sets of rG4detector
data_sets = ["test"]
for data in data_sets:
    data_df = pd.read_csv(f'{data}_data.csv')

    # get position locations in the most prominent transcript
    tr_data_df = pd.DataFrame(columns=["transcript", "position", "rsr", "total_reads"])
    for chrom in data_df["chromosome"].unique():
        for strand in ("+", "-"):
            print(f"Extracting locations for {chrom}{strand}")
            s_data = data_df[(data_df["chromosome"] == chrom) & (data_df["strand"] == strand)]
            s_locations_df = get_transcript_and_position(pos_list=s_data["position"].tolist(),
                                                         tr_list=chr_dict[chrom][strand],
                                                         rsr_list=s_data["rsr"].tolist(),
                                                         reads_list=s_data["total_reads"].tolist(),
                                                         strand=strand)
            tr_data_df = pd.concat([tr_data_df, s_locations_df], ignore_index=True, axis=0)
    tr_data_df.to_csv(f"{data}_tr_data.csv", index=False)

    seqs = []
    f=0
    for idx, row in tr_data_df.iterrows():
        seq_start = row["position"] - SEQ_SIZE - FLANK_SIZE - 1
        seq_end = seq_start + FLANK_SIZE*2 + SEQ_SIZE
        trans_seq = str(transcripts_dict[row["transcript"]].seq)
        seq = trans_seq[max(0, seq_start):min(len(trans_seq), seq_end)]
        seq = "Z"*max(-seq_start, 0) + seq.upper() + "Z"*max(0, seq_end - len(trans_seq))
        assert len(seq) == FLANK_SIZE*2 + SEQ_SIZE, f"ERROR: seq length == {len(seq)}"
        seqs.append(seq)

    # write sequences to file
    with open(output_dir + f"{data}-seq", 'w') as f:
        for s in seqs:
            f.write(f"{s}\n")
    # write also in fasta format
    with open(output_dir + f"{data}.fa", 'w') as f:
        for i, s in enumerate(seqs):
            header = tr_data_df.iloc[i]["transcript"] + "-" + str(tr_data_df.iloc[i]["position"])
            f.write(f">{header}\n{s}\n")
print("All done!")
