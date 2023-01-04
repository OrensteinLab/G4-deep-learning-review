#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:28:56 2020

@author: aliceaubert1
"""
from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)
import sys
import argparse
import numpy as np
import pandas as pd
from keras.utils import to_categorical
from keras.models import load_model


def oneHot(string, win=124):
    if len(string) > win:
        pwin = round(win/2)
        c = round(len(string)/2)
        string = string[c-pwin:c+pwin]
    elif len(string) < win:
        z1 = round((win - len(string))/2)
        z2 = win - z1
        string = 'N'*z1 + string + 'N'*z2
		
    trantab = str.maketrans('ACGTNZ', '012344')
    string = string + 'ACGTN'
    data = list(string.translate(trantab))
    mat = to_categorical(data)[0:-5]
    mat = np.delete(mat, -1, -1)
    mat = mat.astype(np.uint8)
    return mat
  
def read_files(path, window=125):
    f = pd.read_csv(path, header=None, names=['sequence'], index_col=False)
    f[f['sequence'].str.contains('>')] = '>'
    f = pd.DataFrame(''.join(f['sequence'].values).split('>')[1:]).reset_index(drop=True)
    return f


def user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', help='Path to fasta file. of rG4 data', type=str)
    parser.add_argument('-m', '--model', help='Path to model. from G4detector project', type=str)
    parser.add_argument('-o', '--out', help='Output path', type=str)

    args = parser.parse_args()
    arguments = vars(args)

    return arguments

def main():
    args = user_input()
    data = read_files(args['data'])
    x = np.array(list(map(oneHot, list(data[0]))))
    model = load_model(args['model'])
    pred = model.predict(x)
    data['predicted'] = pred.squeeze()
    data['predicted'].to_csv(args['out'], index=False)


if __name__ == "__main__":
    main()
