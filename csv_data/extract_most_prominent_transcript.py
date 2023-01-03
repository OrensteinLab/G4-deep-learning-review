import re
import gzip
import pickle
import os
from time import time
import sys
import pandas as pd

chr_dict = {}


class transcript:
    def __init__(self, name, length, start, end, tsl, gene, strand):
        self.name = name
        self.len = length
        self.exons_ranges = []
        self.start = int(start)
        self.end = int(end)
        self.gene = gene
        self.tsl = 6 if tsl == "NA" else int(tsl)
        self.strand = strand


    def add_exon(self, start, end):
        pos = 0
        start, end = int(start), int(end)
        while pos < len(self.exons_ranges):
            if self.strand == "+":
                if start < self.exons_ranges[pos][0]:
                    break
            else:
                if start > self.exons_ranges[pos][0]:
                    break
            pos += 1
        self.exons_ranges.insert(pos, (start, end))

    def __str__(self):
        print_str = f"name = {self.name}\nlength = {self.len}\nstart = {self.start}\nend = {self.end}\n" \
                    f"tsl = {self.tsl}\nexons = {self.exons_ranges}"
        return print_str


def check_convert_chr_id(chr_id):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.
    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...
    # >>> chr_id = "chrX"
    # >>> check_convert_chr_id(chr_id)
    'chrX'
    # >>> chr_id = "4"
    # >>> check_convert_chr_id(chr_id)
    'chr4'
    # >>> chr_id = "MT"
    # >>> check_convert_chr_id(chr_id)
    'chrM'
    # >>> chr_id = "GL000009.2"
    # >>> check_convert_chr_id(chr_id)
    False
    # >>> chr_id = "chrUn_KI270442v1"
    # >>> check_convert_chr_id(chr_id)
    False
    """
    assert chr_id, "given chr_id empty"

    if re.search("^chr", chr_id):
        if not re.search("^chr[\dMXY]+$", chr_id):
            chr_id = False
    else:
        # Convert to "chr" IDs.
        if chr_id == "MT":
            chr_id = "M"
        if re.search("^[\dMXY]+$", chr_id):
            chr_id = "chr" + chr_id
        else:
            chr_id = False
    return chr_id


def gff_get_transcript_lengths(in_gff, tr2exc_dic=None):
    """
    Get transcript lengths (= length of their exons, not unspliced length!)
    from GFF file.
    tr2exc_dic:
    Optionally provide a transcript ID to exon count dictionary for counting
    transcript exons.
    # >>> in_gtf = "test_data/map_test_in.gtf"
    # >>> gtf_get_transcript_lengths(in_gtf)
    {'ENST001': 2000, 'ENST002': 2000}
    """

    t = time()
    # Transcript ID to exonic length dictionary.
    tr2len_dic = {}
    # Open GFF either as .gz or as text file.
    if re.search(".+\.gz$", in_gff):
        fp = gzip.open(in_gff, 'rt')
    else:
        fp = open(in_gff, "r")
    for line in fp:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        infos = cols[8]
        if not feature == "exon":
            continue
        # Extract transcript ID.
        m = re.search('transcript_id=(.+?);', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" % (in_gff, line)
        tr_id = m.group(1)
        # Sum up length.
        ex_len = feat_e - feat_s + 1
        if tr_id not in tr2len_dic:
            tr2len_dic[tr_id] = ex_len
        else:
            tr2len_dic[tr_id] += ex_len
        if tr2exc_dic is not None:
            if tr_id not in tr2exc_dic:
                tr2exc_dic[tr_id] = 1
            else:
                tr2exc_dic[tr_id] += 1
    fp.close()
    assert tr2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_gff)

    print(f"gff_get_transcript_lengths execution time = {round(time()-t)}s")
    return tr2len_dic


def sort_transcripts():
    for chrom in chr_dict:
        for strand in chr_dict[chrom]:
            chr_list = []
            for tr in chr_dict[chrom][strand]:
                pos = 1
                while pos <= len(chr_list):
                    if chr_dict[chrom][strand][tr].start >= chr_list[-pos].start:
                        break
                    pos += 1
                chr_list.insert(len(chr_list)-pos+1, chr_dict[chrom][strand][tr])
            chr_dict[chrom][strand] = chr_list


def get_transcripts_locations(in_gff, out_dir):
    id2sc = {}
    for i in range(5):
        pos = i + 1
        pos_str = "%i" % pos
        id2sc[pos_str] = pos
    id2sc["NA"] = 6

    # Read in transcript length (exonic regions).
    print("Read in transcript lengths (exonic lengths) from GTF ... ")
    tr2exc_dic = {}
    tr2len_dic = gff_get_transcript_lengths(in_gff, tr2exc_dic=tr2exc_dic)
    assert tr2len_dic, "no transcript lengths read in from --gtf (invalid file format?)"
    print("# transcripts read in:  %i" % (len(tr2len_dic)))



    print("Extract most prominent transcripts ... ")

    # Open GFF either as .gz or as text file.
    if re.search(".+\.gz$", in_gff):
        f = gzip.open(in_gff, 'rt')
    else:
        f = open(in_gff, "r")

    t = time()
    for i, line in enumerate(f):
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        # Take only transcripts and exons
        if feature not in ("transcript", "exon"):
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id)
        if not new_chr_id:
            continue

        m = re.search('transcript_id=(.+?);', infos)
        assert m, "transcript_id entry missing in GFF file \"%s\", line \"%s\"" % (in_gff, line)
        tr_id = m.group(1)
        m = re.search('gene_id=(.+?);', infos)
        assert m, "transcript_id entry missing in GFF file \"%s\", line \"%s\"" % (in_gff, line)
        gene_id = m.group(1)

        # Transcript length.
        tr_len = tr2len_dic[tr_id]

        # Get transcript support level (TSL).
        m = re.search('transcript_support_level=(.+?);', infos)
        tsl_id = "NA"
        if m:
            tsl_id = m.group(1)

        # Add transcript
        if feature == "transcript":
            if chr_id not in chr_dict:
                chr_dict[chr_id] = {"+": {}, "-": {}}
            chr_dict[chr_id][feat_pol][tr_id] = transcript(tr_id, tr_len, feat_s, feat_e, tsl_id, gene_id, feat_pol)
        else:
            chr_dict[chr_id][feat_pol][tr_id].add_exon(feat_s, feat_e)
    f.close()

    print(f"Done reading lines - {round((time()-t))}s")

    print("Sorting transcripts")
    t = time()
    sort_transcripts()
    print(f"Done Sorting - {round((time()-t))}s")


    with open(f'{out_dir}/transcripts_locations.pkl', 'wb') as f:
        pickle.dump(chr_dict, f)


def get_transcript_and_position(pos_list, tr_list, rsr_list, reads_list, strand):
    """
     Assumption: chromosome coordinates are sorted in ascending order
    """
    prominent_tr_idx = None
    tr_idx = 0
    pos_idx = 0
    no_match_sum = 0
    rts_locations = []

    first_match_idx = None
    while pos_idx < len(pos_list):
        # stop if transcript list ended
        if tr_idx >= len(tr_list):
            no_match_sum += len(pos_list[pos_idx:])
            break

        # If current transcript starts after current position
        if tr_list[tr_idx].start > pos_list[pos_idx]:
            # If current position had already match transcript
            if first_match_idx is not None:
                start_point = 0
                found_exon = False
                for exon_range in tr_list[prominent_tr_idx].exons_ranges:
                    if not exon_range[0] < pos_list[pos_idx] < exon_range[1]:
                        start_point += exon_range[1] - exon_range[0] + 1
                    else:
                        found_exon = True
                        if strand == "+":
                            start_point += pos_list[pos_idx] - exon_range[0] + 1
                        else:
                            start_point += exon_range[1] - pos_list[pos_idx] + 1
                        break
                if not found_exon:
                    print("ERROR: Didn't found exon")
                else:
                    spliced = 1 if len(tr_list[prominent_tr_idx].exons_ranges) > 1 else 0
                    rts_locations.append([tr_list[prominent_tr_idx].name, start_point,
                                          rsr_list[pos_idx], reads_list[pos_idx], spliced])
                tr_idx = first_match_idx
                first_match_idx = None
                prominent_tr_idx = None
            # If current position had no matched transcripts
            else:
                no_match_sum += 1
                # print(f"No match for {chrom}-{pos_list[pos_idx]}")
            # next position
            pos_idx += 1
            continue

        # If current transcript ends before current position
        if tr_list[tr_idx].end < pos_list[pos_idx]:
            tr_idx += 1
            continue


        # check if pos is overlapping any exon
        over_lapping_exon = False
        for exon_range in tr_list[tr_idx].exons_ranges:
            if exon_range[0] < pos_list[pos_idx] < exon_range[1]:
                over_lapping_exon = True
                # print("find overlapping trans:")
                # print(tr_list[tr_idx])
                break
        if not over_lapping_exon:
            tr_idx += 1
            continue

        # check for prominent transcript
        if first_match_idx is None:
            first_match_idx = tr_idx
            prominent_tr_idx = tr_idx
        # prominent transcript if tsl is lower, in case of equality, longer is prominent
        elif tr_list[tr_idx].tsl <= tr_list[prominent_tr_idx].tsl and not \
                (tr_list[tr_idx].tsl == tr_list[prominent_tr_idx].tsl and
                 tr_list[tr_idx].len < tr_list[prominent_tr_idx].len):
            prominent_tr_idx = tr_idx

        tr_idx += 1
    print(f"No matched positions = {no_match_sum}/{len(pos_list)}")
    return pd.DataFrame(rts_locations, columns=["transcript", "position", "rsr", "total_reads", "splice"])



if __name__ == "__main__":

    # if len(sys.argv) != 3:
    #     print("ERROR:\nExecution: extract_most_prominent_transcript.py <gff3_path> <output_dir>")
    #     exit(1)
    src = '/data/gencode.v40.primary_assembly.annotation.gff3.gz'
    dst = '/run_mismatch/out'

    # Set output directory
    if not os.path.exists(dst):
        os.makedirs(dst)

    get_transcripts_locations(src, dst)
