# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 16:12:14 2017

@author: Yunjoo Park

G14554-01 Representations and Algorithms for Computational Molecular Biology

Team Project

Ewha Womans University 
Computer Sceince and Engineering
Team G. 
"""

import re
import math
from itertools import izip

states = ('E', 'I')
symbols = ('a', 't', 'g', 'c', 'n')

# get the number of lines in lineCounts
def getCounts(types, filename):
    f = open(types + '/' + filename, 'r')
    counts = [x for x in f.read().splitlines()]
    return counts

# extract [chrom, start, end] from a given BED file
def extract3ElementsFromBED(types, bedFileName):
    print 'Reading from ' + bedFileName + '...'
    bedFile = open(types + '/' + bedFileName + '.csv','r')
    bedByLines = (x for x in bedFile.read().splitlines())
    bedFile.close()
    
    bedElementsByLines = (x[0:3] for x in (line.strip().split(',') for line in bedByLines))
    bed4Elements = ([chrom, int(start), int(end)] for chrom, start, end in bedElementsByLines)
    
    return bed4Elements

# extract [chrom, start, end] from a given BED file
def extract3ElementsFromBEDwithTab(types, bedFileName):
    #print 'Reading from ' + bedFileName + '...'
    bedFile = open(types + '/' + bedFileName + '.csv','r')
    bedByLines = (x for x in bedFile.read().splitlines())
    bedFile.close()
    bedElementsByLines = (x[0:3] for x in (line.strip().split('\t') for line in bedByLines))
    bed4Elements = ([chrom, int(start), int(end)] for chrom, start, end in bedElementsByLines)
    return bed4Elements


# return chromosome whole sequence from fasta
def chr_whole_seq(chrom_num):
    f = open('hg19/' + chrom_num + '.fa', 'r')
    #chrString = f.read().lower()
    chr_string = f.read()
    #delete header
    chr_string = re.sub(r'^>.{2,40}','',chr_string)
    #delete white space
    chr_string = re.sub(r'[^natgcNATGC]','',chr_string)
    f.close()
    return chr_string

# return initializing frequency dictionary composed of 'atgcn'
def make_freq_dic():
    base = ['a','t','g','c','n']
    freq_dic = {}
    for x in base:
        freq_dic[x] = 0
    return freq_dic

# build emission dictionary composed of 'atgcn' (sum = 1)
def build_emit_dic(frequency_dic):
    base = ['a','t','g','c','n']
    emit_dic = {}
    total = sum(frequency_dic.values())
    for x in base:
        emit_dic[x] = float(frequency_dic[x]/float(total))
    return emit_dic

# return emission probabilities
# @param data type (exon / intron)
def get_emit_prob(types):
    freq_dic = make_freq_dic()

    linecounts = getCounts(types, 'lineCounts.txt')
    for i in range(len(linecounts)-1):
        
        chromosome = 'chr'+str(i+1)
        elements = extract3ElementsFromBED(types, chromosome)
        element = elements.next()
        whole_seq = chr_whole_seq(element[0])

        #print obs_seq
        print "training...",
        
        whole_seq = whole_seq.lower()
        line_num = 1
        
        for element in elements:
            obs_seq = whole_seq[element[1]-1:element[2]]
            #chr_seq_cpg_len += len(chr_seq_cpg)
            for base in freq_dic:
                #print base,
                count = freq_dic.get(base) + len(re.findall(base, obs_seq))
                #print len(re.findall(base,obs_seq))
                freq_dic.update({base:count})
            line_num += 1
  
        print line_num       
        
    emit_prob = build_emit_dic(freq_dic)
    
    return emit_prob, freq_dic
    
# decoding for a given seqeuence
def viterbi(obs_seq, states, start_prob, trans_prob, emis_prob):
    obs_seq_len = len(obs_seq)
    if obs_seq_len == 0:
        print 'length is 0'
        return 
    V = [{}]
    for state in states:
        V[0][state] = {"prob": math.log(start_prob[state]) + math.log( emis_prob[state][obs_seq[0]]), "prev":None}

    for t in range(1, len(obs_seq)):
         V.append({})
         for cur_st in states:
             max_tr_prob = max(V[t-1][prev_st]["prob"] + math.log(trans_prob[prev_st][cur_st]) for prev_st in states)
             for prev_st in states:
                 if (V[t-1][prev_st]["prob"] + math.log( trans_prob[prev_st][cur_st])) == max_tr_prob:
                     max_prob = max_tr_prob + math.log( emis_prob[cur_st][obs_seq[t]])
                     V[t][cur_st] = {"prob": max_prob, "prev": prev_st}
                     break
    """
    for line in dptable(V):
        print line
    """
    # Calculate the highest probability
    max_prob = max(value["prob"] for value in V[-1].values())
    #print V
    
    previous = None
       
    opt = []
    # Get the most probable state and its backtrack
    for state, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(state)
            previous = state
            break
    
    # backtracking from the last to the first            
    for i in range(len(V)-2, -1, -1):
        opt.insert(0, V[i+1][previous]["prev"])
        previous = V[i+1][previous]["prev"]
        
    print 'The steps of states are'
    print ''.join(opt)
    print 'with the hightest probabilty of %s' % max_prob

# Print a table of steps from dictionary
def dptable(V):
     yield " ".join(("%12d" % i) for i in range(len(V)))
     for state in V[0]:
         yield "%.7s: " % state + " ".join("%.7s" % ("%f" % v[state]["prob"]) for v in V)

# Return the start probability of the given state
# if 'state' is not in the state set of the model, return 0
def start_prob_fun(state):
    if state not in states:
        print 'state not in states start_prob_fun',
        return 0
    return start_prob[state]

# Return the probability that the transition from 'state_from' to 'state_to'
# If either 'state_from' or 'state_to' are not contained in the staet set of
# the model, return 0
def trans_prob_fun(state_from, state_to):
    if state_from not in states or state_to not in states:
        print 'state not in states trans_prob_fun',
        return 0
    return trans_prob[state_from][state_to]

# Return the emission probability for 'symbol' associated with the 'state'
# If either 'state' or 'symbol' are not contained in the model, return 0
def emit_prob_fun(state, symbol):
    if state not in states or symbol not in symbols:
        print 'state or symbol not in states/symbols emit_prob_fun',
        return 0
    return emit_prob[state][symbol]

# forward algorithm
# @return probability
def forward(sequence):
    sequence_length = len(sequence)
    if sequence_length == 0:
        return []

    alpha = [{}]
    for state in states:
        alpha[0][state] = start_prob[state] * emit_prob[state][sequence[0]]
        
    for index in xrange(1, sequence_length):
        alpha.append({})
        for state_to in states:
            prob = 0
            for state_from in states:
                prob += alpha[index - 1][state_from] * trans_prob[state_from][state_to]
            alpha[index][state_to] = prob * emit_prob[state_to][sequence[index]]

    return alpha

# backward 
def backward(sequence):
    sequence_length = len(sequence)
    if sequence_length == 0:
        return []

    beta = [{}]
    for state in states:
        beta[0][state] = 1

    for index in xrange(sequence_length - 1, 0, -1):
        beta.insert(0, {})
        for state_from in states:
            prob = 0
            for state_to in states:
                prob += beta[1][state_to] * trans_prob_fun(state_from, state_to) * emit_prob_fun(state_to, sequence[index])
            beta[0][state_from] = prob

    return beta

# Return probability for a given sequence
# Use the 'forward algorithm' to evaluate the given sequence.
def evaluate(sequence):
    length = len(sequence)
    
    # return 0 when sequence length = 0
    if length == 0:
        return 0

    prob = 0.0
    
    # get alpha using 'forward' algorithm
    alpha = forward(sequence)
    
    for state in alpha[length - 1]:
        prob += alpha[length - 1][state]

    return prob

# Train a HMM model using the given sequences
# Learning algorithm will stop when the diff of the log-likelihood 
# between two consdecutive iterations is less than delta.
# @param sequeces, delta=0.0001, smoothing=0
# update each probabilities (emission prob, transmission prob)
def train(sequences, delta=0.0001, smoothing=0):
    length = len(sequences) # 데이터 묶음 개수 확인 

    
    old_likelihood = 0 
    for _, symbol_list in sequences:

        old_likelihood += math.log(evaluate(symbol_list)) 

    old_likelihood /= length
    
    #print '\n new likelihood train start \n'
    while True:
        new_likelihood = 0
        
        for _, symbol_list in sequences:
            learn(symbol_list, smoothing) # 한 단위의 시퀀스 학습 
            # check probability
            check_prob = evaluate(symbol_list)

            if check_prob == 0:
                print 'here',
                continue
            else:
##                print 'check_prob'
##                print check_prob
                new_likelihood += math.log(check_prob)

        new_likelihood /= length

        if abs(new_likelihood - old_likelihood) < delta:
            break

        old_likelihood = new_likelihood #교체

# Use the given 'sequence' to find the best state transition and emission probabilities.
# @param sequence, smoothing = 0
# update each probabilities (emission prob, transmission prob)
def learn(sequence, smoothing=0):
    length = len(sequence)
    alpha = forward(sequence)
    beta = backward(sequence)

    gamma = [] #
    for index in range(length):
        prob_sum = 0
        gamma.append({})
        for state in states:
            prob = alpha[index][state] * beta[index][state]
            gamma[index][state] = prob
            prob_sum += prob

        if prob_sum == 0:
            continue

        for state in states:
            gamma[index][state] /= prob_sum

    xi = []
    for index in range(length - 1):
        prob_sum = 0
        xi.append({})
        for state_from in states:
            xi[index][state_from] = {}
            for state_to in states:
                prob = alpha[index][state_from] * beta[index + 1][state_to] * trans_prob_fun(state_from, state_to) * emit_prob_fun(state_to, sequence[index + 1])
                xi[index][state_from][state_to] = prob
                prob_sum += prob

        if prob_sum == 0:
            continue

        for state_from in states:
            for state_to in states:
                xi[index][state_from][state_to] /= prob_sum

    states_number = len(states)
    symbols_number = len(symbols)
    for state in states:
        # update start probability
        #start_prob[state] = (smoothing + gamma[0][state]) / (1 + states_number * smoothing)

        # update transition probability
        gamma_sum = 0
        for index in xrange(length - 1):
            gamma_sum += gamma[index][state]

        if gamma_sum > 0:
            denominator = gamma_sum + states_number * smoothing
            for state_to in states:
                xi_sum = 0
                for index in xrange(length - 1):
                    xi_sum += xi[index][state][state_to]
                    trans_prob[state][state_to] = (smoothing + xi_sum) / denominator
        else:
            for state_to in states:
                trans_prob[state][state_to] = 0
        

        # update emission probability
        gamma_sum += gamma[length - 1][state]

        
        emit_gamma_sum = {}
        for symbol in symbols:
            emit_gamma_sum[symbol] = 0

        for index in xrange(length):
            
            emit_gamma_sum[sequence[index]] += gamma[index][state]## emit_gamma_sum = {}    gamma = []

        if gamma_sum > 0:
            denominator = gamma_sum + symbols_number * smoothing
            for symbol in symbols:
                emit_prob[state][symbol] = (smoothing + emit_gamma_sum[symbol]) / denominator
        else:

            for symbol in symbols:
                emit_prob[state][symbol] = 0

def init():
    print '> Initialize probabilities'
    trans_prob = {'E': {'E':0.9, 'I':0.1}, 'I': {'I':0.9, 'E':0.1}}
    start_prob = {'E':0.8, 'I':0.2}
    emit_prob_exon, freq_dic_exon = get_emit_prob('exon')
    emit_prob_intron, freq_dic_intron = get_emit_prob('intron')
    emit_prob = {'E':emit_prob_exon, 'I':emit_prob_intron}
    """
    print '1. start probaility'
    print start_prob
    print '2. transition probability'
    print trans_prob
    print '3. emission probability'
    print emit_prob
    """
    return start_prob, trans_prob, emit_prob, freq_dic_exon, freq_dic_intron


def make_trainingdata(i):
    print '> Make training data'
    state_list = []

    chromosome = 'chr'+str(i) 
    exon_elements = extract3ElementsFromBED('exon', chromosome)
    counts = getCounts('entire', chromosome + '_exon_counts.txt')
    whole_seq = chr_whole_seq(chromosome).lower()
    linecounts = getCounts('entire', 'lineCounts.txt')
    #print linecounts[i-1]
    for i in range(int(linecounts[i-1])):
        state_ = []
        intron_start = 0
        intron_end = 0
        intron_length = 0
        exon_element = exon_elements.next()
        exon_length = exon_element[2] - exon_element[1] + 1
        intron_start = exon_element[2] + 1
        #print 'exon length : ' + str(exon_element[1]) + ' , '  + str(exon_element[2])
        state_.extend( ['E' for x in range(0, exon_length)])
        for j in range(1, int(counts[i])):
                exon_element = exon_elements.next()
                exon_length = exon_element[2] - exon_element[1] + 1
                
                intron_end = exon_element[1] - 1
                intron_length = (intron_end)  - (intron_start) + 1
                state_.extend( ['I' for x in range(0, intron_length)])
                
                #print 'intron_length : ' + str(intron_start) + ' , ' + str(intron_end)
                #print 'exon_length : ' + str(exon_element[1]) + ' , ' + str(exon_element[2])
                intron_start = exon_element[2] + 1

                state_.extend( ['E' for x in range(0, exon_length)])
        
        state_list.append(state_)
        
        # make symbol list
        elements = extract3ElementsFromBEDwithTab('entire', 'chr22')
        obs_seq_list = []
        for i in range(830):
            element = elements.next()
            obs_list = list(whole_seq[element[1]-1:element[1]+299])
            obs_seq_list.append(obs_list)
            #length = len(obs_seq)
            
            #print str(len(training_state_list[i][0:300])) + ' seq len: ' + str( length)
            
        training_data = []
        for state_item, symbol_item in izip(state_list, obs_seq_list):
            combine =[]
            combine.append(state_item)
            combine.append(symbol_item)
            training_data.append(combine)

        return training_data




###################################
# Initializing probabilites
###################################
start_prob, trans_prob, emit_prob, freq_dic_exon, freq_dic_intron = init()
print '> Results'
###################################
# testing before learning (with initial probabilities)
###################################
print 'before----\ntrans_prob'
print trans_prob
print 'emit_prob'
print emit_prob
elements = extract3ElementsFromBEDwithTab('entire', 'chr22')
# read whole sequence in fasta
whole_seq = chr_whole_seq('chr22').lower()
n = 1

for i in range(10):
    element = elements.next()

    obs_seq = whole_seq[element[1]-1:element[1]+299]
    
    length = len(obs_seq)
    print 'Line ' + str(n) + ': chr22 from ' + str(element[1]) + ' to ' + str(element[1]+300)
    
    print '>>> observed sequence:',
    print obs_seq[0:10] + '...' + obs_seq[-10:-1] + ' with length ' + str(length)
    print '>>> Evaluation: %e' % evaluate(list(obs_seq))
    print '>>> Decoding: '
    viterbi(obs_seq, states, start_prob, trans_prob, emit_prob)
    n += 1


###################################
# training 
###################################
#symbols = ('a', 't', 'g', 'c', 'n')
print '\n>training..........'
#for i in range(1,22):
training_data = make_trainingdata(21)
train(training_data)
    
###################################
# testing after learning
###################################
print '\nafter----\ntrans_prob'
print trans_prob
print 'emit_prob'
print emit_prob

elements = extract3ElementsFromBEDwithTab('entire', 'chr22')
# read whole sequence in fasta
whole_seq = chr_whole_seq('chr22').lower()
n = 1

for i in range(10):
    element = elements.next()

    obs_seq = whole_seq[element[1]-1:element[1]+299]
    
    length = len(obs_seq)
    print 'Line ' + str(n) + ': chr22 from ' + str(element[1]) + ' to ' + str(element[1]+300)
    
    print '>>> observed sequence:',
    print obs_seq[0:10] + '...' + obs_seq[-10:-1] + ' with length ' + str(length)
    print '>>> Evaluation: %e' % evaluate(list(obs_seq))

