import numpy as np
from pathlib import Path

def combinations_builder(chars, m):
    n = len(chars)
    indices = [0]*m
    # strings = []
    while True:
        string = "".join(chars[i] for i in indices)
        yield string
        for i in range(m-1, -1, -1):
            if indices[i] < n-1 :
                indices[i] += 1
                break
            else :
                indices[i] = 0
        else:
            return

def populate_dict(dict, sequence, m, pseudocount):
    window_size = m+1
    for i in range(len(sequence)- window_size):
        window_seq = sequence[i:i+window_size]
        if window_seq in dict :
            dict[window_seq] += 1
    for key in dict:
        dict[key] += pseudocount
    return dict

def markov_model(string_freq1,string_freq2, chars):
    transition_matrix = np.zeros((len(string_freq2), len(chars)))

    from_seq_list = list(string_freq2.keys())

    for i in range(len(string_freq2)):
        from_seq = from_seq_list[i]
        for j in range(len(chars)):
            complete_seq = from_seq + chars[j]
            count = string_freq1[complete_seq]
            nf = string_freq2[from_seq]
            transition_matrix[i, j] = count/nf
    
    return transition_matrix

# print(from_state_idx)
def log_likelihood(transition_matrix, string_freq2, chars, test_seq, m):
    from_state_idx = dict()
    for i, k in enumerate(string_freq2):
        from_state_idx[k] = i
    to_state_idx = dict()
    for i, k in enumerate(chars):
        to_state_idx[k] = i
    window_size = m+1
    log_likelihood_score = 0
    for i in range(len(test_seq)-m):
        window_seq = test_seq[i:i+window_size]
        from_state = window_seq[:-1]
        # print(from_state)
        to_state = window_seq[-1]
        if from_state in from_state_idx and to_state in to_state_idx:
            from_idx = int(from_state_idx[from_state])
            to_idx = int(to_state_idx[to_state])
            transition_prob = transition_matrix[from_idx, to_idx]
            log_likelihood_score += np.log(transition_prob)
    return log_likelihood_score
