import argparse
import csv

#function to load a FASTA file into a dictionary, >ID is used as key, sequence is used as value
def fasta_to_dict(fasta_file):
    fasta_dic = {}
    name = None
    lines = open(fasta_file).readlines()
    for line in lines:
        if line.startswith('>'):
            if name is not None:
                fasta_dic[name] = new_seq
            new_seq = ''
            name = line[1:].strip()
        else:
            seq = line.strip().rstrip()
            new_seq+=seq
    fasta_dic[name] = new_seq
    #print(fasta_dic)
    return fasta_dic

# function to parse sequence strings from txt file into a list, 1 sequence per line
def read_in_querry_sequences(filename):
    sequence_list = []
    lines = open(filename).readlines()
    for line in lines:
        seq = line.strip().rstrip()
        sequence_list.append(seq)
    return sequence_list

#function to find sequences from txt file and find index of sequence in a FASTA file
def sequence_finder(fasta_file, sequence):
    results_dict = {}
    fasta_dict = fasta_to_dict(fasta_file)
    sequence_list = read_in_querry_sequences(sequence)
    for seq in sequence_list:
        for key, value in fasta_dict.items():
            if seq in value:
                if key in results_dict:
                    results_dict[key].append((seq, KMPSearch(seq,value)))
                else:
                    #results_dict[key] = list([seq, KMPSearch(seq,value)])
                    results_dict[key] = [(seq, KMPSearch(seq,value))]
    filednames=['FASTA_ID', 'Sequence_found_and_start_ index']
    with open(fasta_file.split('.')[0] + 'sequence_found.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(filednames)
        for k, v in results_dict.items():
            writer.writerow([k, *v]) # *expands list in value
    return results_dict



def find_all_substrings(a_str, sub):
    start =0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)

def KMPSearch(pat, txt):
    index = []
    M = len(pat)
    N = len(txt)
    # lps array must be as long as pattern
    lps = [0]*M
    j = 0 # index for pattern
    # calculate lps[] array used in KMPSearch
    computeLPSArray(pat, M, lps)

    i = 0 # index for text to be searched
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1
        if j == M:
            index.append(i-j)
            #index += ' '
            j = lps[j-1]

        # mismatch after j matches
        elif i < N and pat[j] != txt[i]:
            # No need to back-track more than one character is lps[]
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return index
#function to build a longest prefix suffix array needed for KMP search to work
def computeLPSArray(pat, M, lps):
    len = 0 # 1st pointer in pattern = length of the previous longest prefix suffix

    lps[0] =0 # one character cannot have a suffix, so lps[0] is always 0
    i = 1 #2nd pointer in pattern

    # loop over length of pattern
    while i < M:
        if pat[i]== pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            if len != 0:
                len = lps[len-1]

            # When pat[i] !=pat[len] len is not incremented.
            else:
                lps[i] = 0
                i += 1

#argumnet parsing and help
parser = argparse.ArgumentParser(description='This script finds the index of a list of substrings in a FASTA file and returns the FASTA header, number of instances and start index for each match')
parser.add_argument('fasta_path', type=str, help='absolute path to FASTA file to be searched')
parser.add_argument('seq_path', type=str, help='absolute path to text file with search sequences')
args = parser.parse_args()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    path = args.seq_path
    fasta_path = args.fasta_path
    # print(read_in_querry_sequences(path))
    print(sequence_finder(fasta_path,path))

