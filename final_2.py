import numpy as np
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import argparse 
import tqdm as tqdm
import matplotlib.pyplot as plt
import sklearn.metrics as met
from functions import *
import time
import tracemalloc
import pandas as pd
from datetime import datetime

np.random.seed(42)

def the_ultimate_function(tsv_path, tf_name, m = 1, k = 3):
    tsv_file = Path(tsv_path)

    TF_indexes = [3,4,5,6]
    
    TF_names = []
    with open(tsv_file, "r") as f:
        # next(f, None)
        binding_info = []
        for i,l in enumerate(f):
            l = l.strip()
            elements = l.split('\t')
            if i == 0:
                TF_names.extend([elements[j] for j in TF_indexes])

    TF_dict = dict(zip(TF_names,TF_indexes))

    TF_idx = TF_dict[tf_name]

    with open(tsv_file, "r") as f:
        next(f, None)
        binding_info = []
        for l in f:
            elements = l.split('\t')
            u_b_string = "".join(elements[TF_idx])
            binding_info.append(u_b_string.strip())

    binding_info = np.array(binding_info)
    binding_info_dict = {'U': np.where(binding_info == 'U')[0], 
                        'B': np.where(binding_info == 'B')[0]}

    # print(binding_info_dict)

    """
    Reading Fasta
    """

    fasta_file = Path(f"./fasta_outputs/{tsv_file.stem}.fa")

    seq_dict = dict()

    with open(fasta_file , 'r') as ff:
        seq_idx = -1
        for l in ff:
            l = l.strip()
            if l.startswith('>'):
                seq_idx += 1
            else:
                if seq_idx in seq_dict:
                    seq_dict[seq_idx] += l
                else:
                    seq_dict[seq_idx] = l
    # print(seq_dict)
    # k-cross validation

    # m = 1
    chars = ['A', 'T', 'G', 'C']
    # possible_combinations = []
    pseudocount = 1

    # k = 3

    U_idx = binding_info_dict['U']
    B_idx = binding_info_dict['B']
    np.random.shuffle(U_idx)
    np.random.shuffle(B_idx)
    U_split = np.array_split(U_idx, k)
    B_split = np.array_split(B_idx, k)        

    log_dict = {
                'chrm':[],
                'TF':[],
                'order':[],
                'fold_idx':[],
                'start_time':[],
                'stop_time':[],
                'time_elapsed':[],
                'AOC-ROC':[],
                'AOC-PR':[],
                'memory_used (MB)':[]
                }
    for i in tqdm.tqdm(range(k)):
        tracemalloc.start()
        tic = time.perf_counter()
        start_time = datetime.now()
        train_indices = [j for j in range(k) if j != i]
        # print(train_indices)
        train_block_pos = set(np.concatenate([B_split[l] for l in train_indices]))
        train_block_neg = set(np.concatenate([U_split[l] for l in train_indices]))
        test_block_pos = set(B_split[i])
        test_block_neg = set(U_split[i])

        # lets create the dictionary of characters
        string_freq1_pos = {s: 0 for s in combinations_builder(chars, m+1)}
        string_freq2_pos = {s: 0 for s in combinations_builder(chars, m)}
        string_freq1_neg = string_freq1_pos.copy()
        string_freq2_neg = string_freq2_pos.copy()

        # training part
        tqdm.tqdm.write(f"training started for fold-{i}")
        # for seq_idx in tqdm.tqdm():
        
        for seq_idx in tqdm.tqdm(train_block_pos):
            # print(True)
            sequence = seq_dict[seq_idx]
            populate_dict(string_freq1_pos, sequence, m, pseudocount)
            populate_dict(string_freq2_pos, sequence, m-1, pseudocount)
        for seq_idx in tqdm.tqdm(train_block_neg):
            # print(True)
            # print(seq_idx, sequence)
            sequence = seq_dict[seq_idx]
            populate_dict(string_freq1_neg, sequence, m, pseudocount)
            populate_dict(string_freq2_neg, sequence, m-1, pseudocount)        
        tm_pos = markov_model(string_freq1_pos, string_freq2_pos, chars)
        tm_neg = markov_model(string_freq1_neg, string_freq2_neg, chars)
        # print(tm_pos)
        # print(tm_neg)
        tqdm.tqdm.write(f"training completed for fold-{i}...starting test")
        log_likelihood_scores = []
        y_true = []
        # ll_pos_bits = []
        # ll_neg_bits = []
        seq_len = len(sequence)

        # for idx, seq in tqdm.tqdm(seq_dict.items()):
            # sequence_2 = seq_dict[idx]
        for idx in tqdm.tqdm(test_block_pos | test_block_neg):
            seq = seq_dict[idx]
            ll_pos = log_likelihood(tm_pos, string_freq2_pos, chars, seq, m)
            # ll_pos_bits.append(ll_pos/seq_len)
            ll_neg = log_likelihood(tm_neg, string_freq2_neg, chars, seq, m)
            # ll_neg_bits.append(ll_neg/seq_len)
            score = float(ll_pos) - float(ll_neg)
            # print(score)
            log_likelihood_scores.append(score)

            if idx in test_block_pos:
                y_true.append(1)
            else:
                y_true.append(0)

        # print(log_likelihood_scores)

        fpr, tpr, thresholds = met.roc_curve(y_true, log_likelihood_scores)
        roc_auc = met.roc_auc_score(y_true, log_likelihood_scores)
        precision, recall, thresholds = met.precision_recall_curve(y_true, log_likelihood_scores)
        pr_auc = met.average_precision_score(y_true, log_likelihood_scores)
        tqdm.tqdm.write(f"Area Under Curve (AUC): {roc_auc:.4f}")
        tqdm.tqdm.write(f"AUC-PR (Average Precision): {pr_auc:.4f}")
        toc = time.perf_counter()
        end_time = datetime.now()
        t = toc - tic
        print('time taken to calculate auc-roc and auc-pr', t, "in seconds")
        baseline = sum(y_true) / len(y_true)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # ROC curve : Gemini AI model was used here to obtain the code blocks for plotting purposes
        ax1.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
        ax1.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Random Guess (0.5)')
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlabel('False Positive Rate', fontsize=12)
        ax1.set_ylabel('True Positive Rate', fontsize=12)
        ax1.set_title('Receiver Operating Characteristic (ROC)', fontsize=14)
        ax1.legend(loc="lower right")
        ax1.grid(True, alpha=0.3)

        # Precision-Recall Curve
        ax2.plot(recall, precision, color='purple', lw=2, label=f'PR curve (AP = {pr_auc:.3f})')
        ax2.axhline(y=baseline, color='gray', linestyle='--', label=f'Baseline ({baseline:.2f})')
        ax2.set_xlim([0.0, 1.0])
        ax2.set_ylim([0.0, 1.05])
        ax2.set_xlabel('Recall (Sensitivity)', fontsize=12)
        ax2.set_ylabel('Precision', fontsize=12)
        ax2.set_title('Precision-Recall Curve', fontsize=14)
        ax2.legend(loc="upper right")
        ax2.grid(True, alpha=0.3)

        # Save and Close
        plt.tight_layout()
        output_file = Path(f"./plots/metrics_plot_{m}_{tsv_file.stem}_{tf_name}_fold_{i}.png")
        plt.savefig(output_file, dpi=300)
        print(f"Combined plot saved to {output_file}")
        plt.close()
        log_dict['chrm'].append(tsv_file.stem[:4])
        log_dict['order'].append(m)
        log_dict['TF'].append(tf_name)
        log_dict['fold_idx'].append(f"{i}/{k}")
        log_dict['start_time'].append(start_time.time())
        log_dict['stop_time'].append(end_time.time())
        log_dict['time_elapsed'].append(round(t, 3))
        log_dict['AOC-ROC'].append(round(roc_auc, 3))
        log_dict['AOC-PR'].append(round(pr_auc, 3))
        _,peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        module_memory = tracemalloc.get_tracemalloc_memory()
        log_dict['memory_used (MB)'].append((peak-module_memory)/(1024**2))
    return log_dict
    # return roc_auc, pr_auc

def main():
        parser = argparse.ArgumentParser(description="Markov Model Cross Validation")
        parser.add_argument('--order', type = int, default=1, help="Order of the Markov model(m)")
        parser.add_argument('--k', type = int, default=3, help = "k-value")
        parser.add_argument('--input', type = Path, required=True, help = 'enter the path of the tsv file')
        parser.add_argument('--TF', type = str, required=True, default='CTCF', help = 'enter the name of the transcription factor')
        # parse args
        args = parser.parse_args()
        m = args.order
        k = args.k
        tsv_file = Path(args.input)
        tf = args.TF
        print(f"Running with: Order={m}, K={k}, TF={tf}, File={tsv_file.name}")
        a = the_ultimate_function(tsv_file, tf, m, k)
        print(a)
        # df.to_csv(f"/home/photon/CFG-course_project/log_m_{m}_k_{k}_{tf}.txt", sep = '\t', header = True)
# the know-how to write the code block below was learnt from the internet reading various articles and using Gemini to understand some parts
if __name__ == '__main__':
    # try:
    main()
    # except Exception as e:
    #     print("An error occured:", e)
    #     import traceback
    #     traceback.print_exc