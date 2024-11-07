import sys
def read_fasta(file):
    with open(file, 'r') as f:
        sequences = []
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        sequences.append(sequence)
    #print (sequences)
    return sequences

def save_result(output_file, align1, align2, align_score):
    result = f"align seq's:\n{align1}\n{align2}\n\nscore: {align_score}\n"

    with open(output_file, 'w') as file:
        file.write(result)

def needleman_wunsch(seq1, seq2, match, mismatch, gap):
    len_seq1, len_seq2 = len(seq1), len(seq2)
    scores = [[0] * (len_seq2 + 1) for _ in range(len_seq1 + 1)]

    for y in range(1, len_seq1 + 1):
        scores[y][0] = scores[y - 1][0] + gap
    for x in range(1, len_seq2 + 1):
        scores[0][x] = scores[0][x - 1] + gap

    for y in range(1, len_seq1 + 1):
        for x in range(1, len_seq2 + 1):
            match_mismatch = match if seq1[y - 1] == seq2[x - 1] else mismatch
            diagonal = scores[y - 1][x - 1] + match_mismatch
            up = scores[y - 1][x] + gap
            left = scores[y][x - 1] + gap
            scores[y][x] = max(diagonal, up, left)

    align1, align2 = "", ""
    y, x = len_seq1, len_seq2

    while y > 0 and x > 0:
        current_score = scores[y][x]
        if current_score == scores[y - 1][x - 1] + (match if seq1[y - 1] == seq2[x - 1] else mismatch):
            align1 = seq1[y - 1] + align1
            align2 = seq2[x - 1] + align2
            y -= 1
            x -= 1
        elif current_score == scores[y - 1][x] + gap:
            align1 = seq1[y - 1] + align1
            align2 = "-" + align2
            y -= 1
        elif current_score == scores[y][x - 1] + gap:
            align1 = "-" + align1
            align2 = seq2[x - 1] + align2
            x -= 1

    while y > 0:
        align1 = seq1[y - 1] + align1
        align2 = "-" + align2
        y -= 1
    while x > 0:
        align1 = "-" + align1
        align2 = seq2[x - 1] + align2
        x -= 1

    align_score = scores[len_seq1][len_seq2]
    return align1, align2, align_score

def main():
    if len(sys.argv) != 2:
        sys.exit(1)

    ex_fasta = sys.argv[1]
    sequences = read_fasta(ex_fasta)

    match = 1
    mismatch = -1
    gap = -2

    if len(sequences) != 2:
        print("error: FASTA - do not try for less or more than 2seq's")
        sys.exit(1)

    seq1, seq2 = sequences
    align1, align2, align_score = needleman_wunsch(seq1, seq2, match, mismatch, gap)

    save_result("score.txt", align1, align2, align_score)
    print(f"align seq's:\n{align1}\n{align2}\n\nscore: {align_score}\n--> 'score.txt'")

if __name__ == "__main__":
    main()
