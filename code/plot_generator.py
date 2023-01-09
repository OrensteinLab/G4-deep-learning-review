from sklearn.metrics import auc, roc_curve
from scipy.stats import gaussian_kde, pearsonr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

NEGATIVE_TYPES = ['dishuffle', 'pq', 'random']

'''
Figures and Pearson correlation of G4detector (DNA) trained model on rG4 (RNA) test set
'''

def calc_g4detector_on_rg4detector(origin_ds, preds_ds, env, neg):
    origin = pd.read_csv(origin_ds, header=None, names=["transcript", "position", "rsr", "total_reads", "splice"],
                         sep=',')
    preds_g4detector = pd.read_csv(preds_ds, header=None, names=["predicted"], sep=',')
    origin_score = np.array(origin['rsr'][1:].astype(float))
    g4detector_prob = np.array(preds_g4detector['predicted'][1:].astype(float))
    pearson_cor = pearsonr(origin_score, g4detector_prob)
    print(f'{env}_{neg}: pearson: {pearson_cor[0]}, p-value: {pearson_cor[1]}')
    return origin_score, g4detector_prob


def plt_scatter_density(data, origin_score, data_path):
    for type, g4detector_prob in data.items():
        env, neg = type.split('_')
        label = f'{env} {neg}'
        plt.scatter(origin_score, g4detector_prob, label=label, s=10, alpha=0.4)
        xy = np.vstack([origin_score, g4detector_prob])
        z = gaussian_kde(xy)(xy)
        fig, ax = plt.subplots()
        density = ax.scatter(origin_score, g4detector_prob, c=z, s=10, alpha=0.9)
        fig.colorbar(density, label='Density')
        plt.title(label)
        plt.xlabel("Observed RSR ratio")
        plt.ylabel("Predicted G4detector probability")
        pearson_cor = pearsonr(origin_score, g4detector_prob)
        plt.legend([str(pearson_cor[0])[:5]], loc='best', frameon=False, markerscale=0)
        plt.savefig(data_path + f'rsr-g4prob-K-{neg}.jpg')
        plt.show()


def get_g4_scores(env, data_path):
    scores, origin_score = {}, []
    for neg in NEGATIVE_TYPES:
        rg4_data = data_path + 'test_tr_data.csv'
        g4_data = data_path + f'out_g4detector_on_rg4_data/G4detector_prediction_{env}_{neg}.csv'
        origin, pred = calc_g4detector_on_rg4detector(rg4_data, g4_data, env, neg)
        if len(origin_score) == 0:
            origin_score = origin
            scores.update({f'{env}_{neg}': pred})
        else:
            scores.update({f'{env}_{neg}': pred})
    return scores, origin_score


def run_plot_g4detector_on_rg4detector(data_path):
    scores_k, origin_score = get_g4_scores('K', data_path)
    plt_scatter_density(scores_k, origin_score, data_path)

'''
Figures and Pearson correlation of rG4detector (RNA) trained model on G4 (DNA) test set
'''
def run_plot_rg4detector_on_g4detector(data_path):
    scores_k = get_rg4_scores('K', data_path)
    scores_pds = get_rg4_scores('KPDS', data_path)
    plt_auroc(scores_k, 'K', data_path)
    plt_auroc(scores_pds, 'KPDS', data_path)

def get_rg4_scores(env, data_path):
    scores = {}
    for neg in NEGATIVE_TYPES:
        origin_ds = data_path + f'data_out_g4detector_label/{env}/{neg}/test.csv'
        pred_ds = data_path + f'out_rg4detector_on_g4_data/rG4detector_prediction_{env}_{neg}.csv'
        scores.update({f'{env}_{neg}': calc_rg4detector_on_g4detector(origin_ds, pred_ds, env, neg)})
    return scores

def calc_rg4detector_on_g4detector(origin_ds, preds_ds, env, neg):
    origin = pd.read_csv(origin_ds, header=None, names=["seq", "label"], sep=',')
    preds = pd.read_csv(preds_ds, header=None, names=["description", "rG4detector prediction"], sep=',')
    preds_score = preds["rG4detector prediction"][1:].astype(float)
    origin_label = origin["label"][1:].astype(float)
    pearson_cor = pearsonr(origin_label, preds_score)
    print(f'{env}_{neg}: pearson: {pearson_cor[0]}, p-value: {pearson_cor[1]}')
    fp, tp, threshold = roc_curve(origin_label, preds_score)
    au = auc(fp, tp)
    return fp, tp, au

def plt_auroc(data, environment, path):
    for type, values in data.items():
        env, neg = type.split('_')
        fp, tp, au = values
        plt.figure(1)
        plt.plot([0, 1], [0, 1], 'k--')
        l = f'{env} {neg}'
        plt.plot(fp, tp, label= l + '(area = {:.3f})'.format(au))
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
    plt.legend(loc='best')
    plt.savefig(path+f'{environment}_base_auroc.jpg')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-rg4", dest="rg4_path", help="path to g4 data and rg4 predicts")
    parser.add_argument("-g4", dest="g4_path", help="path to rg4 data and g4 predicts")
    args = parser.parse_args()

    run_plot_rg4detector_on_g4detector(args.rg4_path)
    run_plot_g4detector_on_rg4detector(args.g4_path)