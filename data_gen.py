import pandas as pd
import numpy as np


def read_files(pos_f, neg_f):
    # reading from files
    pos = pd.read_csv(pos_f, header=None)[0]
    neg = pd.read_csv(neg_f, header=None)[0]

    # filtering out_g4detector_on_rg4_data chromosome lines
    chrom = pos[pos.str.contains(">")]

    # splitting into chromosome number and window
    chrom = chrom.str.split(':', expand=True)

    header = ['chrom', 'window']
    chrom.columns = header

    # selecting only chromosome number and reseting index
    chrom = chrom.chrom.reset_index(drop=True)

    # filtering out_g4detector_on_rg4_data sequence lines and reseting index
    reads = pos[~pos.str.contains(">")]
    reads = reads.reset_index(drop=True)

    # making new dataframe of positive examples
    pos = pd.DataFrame({'seq': reads, 'chrom': chrom})

    # repeating process for negative examples
    chrom = neg[neg.str.contains(">")]
    chrom = chrom.str.split(':', expand=True)
    header = ['chrom', 'window']
    chrom.columns = header
    chrom = chrom.chrom.reset_index(drop=True)

    reads = neg[~neg.str.contains(">")]
    reads = reads.reset_index(drop=True)

    neg = pd.DataFrame({'seq': reads, 'chrom': chrom})

    # preparing positive set
    p = pos[pos.chrom != '>chr1'].seq  # filtering out_g4detector_on_rg4_data chr1
    p = p.str.upper()
    p = p[~p.str.contains("N")]
    p.reset_index(drop=True, inplace=True)
    lp = np.ones(p.shape[0])

    # chromosome 1 test set
    cp = pos[pos.chrom == '>chr1'].seq  # selecting chr1
    cp = cp.str.upper()
    cp = cp[~cp.str.contains("N")]
    cp.reset_index(drop=True, inplace=True)
    lcp = np.ones(cp.shape[0])

    # negative set
    n = neg[neg.chrom != '>chr1'].seq  # filtering out_g4detector_on_rg4_data chr1
    n = n.str.upper()
    n = n[~n.str.contains("N")]
    n.reset_index(drop=True, inplace=True)
    ln = np.zeros(n.shape[0])

    # chromosome 1 test set (negative)
    cn = neg[neg.chrom == '>chr1'].seq  # selecting chr1
    cn = cn.str.upper()
    cn = cn[~cn.str.contains("N")]
    cn.reset_index(drop=True, inplace=True)
    lcn = np.zeros(cn.shape[0])

    # join and mix
    x = pd.concat([p, n])
    y = np.hstack([lp, ln])
    y = y.astype('int')
    train = pd.DataFrame({'seq': x, 'label': y})
    train = train.sample(frac=1).reset_index(drop=True)

    test_x = pd.concat([cp, cn])
    test_y = np.hstack([lcp, lcn])
    test_y = test_y.astype('int')
    test = pd.DataFrame({'seq': test_x, 'label': test_y})
    test = test.sample(frac=1).reset_index(drop=True)

    return train, test


if __name__ == '__main__':
    environment = ['KPDS', 'K']
    negative_type = ['dishuffle', 'pq', 'random']
    for env in environment:
        for neg in negative_type:
            pos_f = f'/Users/kelimelech/PycharmProjects/transcripts/g4_data_prep/{env}/pos_ex_{env}_125.fa'
            neg_f = f'/Users/kelimelech/PycharmProjects/transcripts/g4_data_prep/{env}/neg_ex_{env}_{neg}_125.fa'
            train, test = read_files(pos_f, neg_f)
            test['seq'] = 'NNN' + test['seq'].astype(str) + 'NNN'
            out_path = f'/Users/kelimelech/PycharmProjects/transcripts/g4_data_prep/demo/'
            with open(out_path + "test.fa", 'w+') as f:
                for i, s in enumerate(test['seq']):
                    header = ' '
                    f.write(f">{header}\n{s}\n")
            test.to_csv(out_path + '/test.csv', index=False, sep=',')
