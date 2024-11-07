**Needleman-Wunsch Sequence Alignment**

This program implements the Needleman-Wunsch algorithm for global sequence alignment in Python. It reads two sequences from a FASTA file, aligns them, and saves the result into a file.

**Requirements**

Python 3.x

**Usage**

Prepare a FASTA file with exactly two sequences:

>seq 1
CCCGCTTTT
>seq 2
GTTCGGG

**Run the script:**

python needleman_wunsch.py sequences.fasta

∨

python3 needleman_wunsch.py sequences.fasta

The aligned sequences and score will be saved in score.txt.

**Output**

The output file score.txt will contain:
 
align seq's:
CCCGCTTTT
-GTTC-GGG

score: -9
