# dna-pattern-matching

This C program performs pattern matching between DNA sequences.

## Description

The pattern matching can be performed using three different algorithms:
- Brute-force
- [Karp-Rabin](https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm)
- [LCSS](https://en.wikipedia.org/wiki/Longest_common_subsequence_problem)

The Brute-force and Karp-Rabin algorithms can be used to find how many occurrences of a smaller DNA sequence are 
found in a larger DNA sequence.

The LCSS algorithm can be used to find the length of the largest common
[subsequence](https://en.wikipedia.org/wiki/Subsequence) between two DNA sequences. It also finds the distance between 
the two DNA sequences.

- Distance is difined as: d(A, B) = 1 - LCSS(A, B) / min(| A |, | B |)

## Usage

Compilation:

```bash
gcc -o dnaPatternMatching dnaPatternMatching.c -lm
```

Run:

```bash
./dnaPatternMatching <algorithm> <dna_sequence_file> <pattern_file or second_dna_sequence_file>
```

Where:

| algorithm |                       |
| --------- | --------------------- |
| -bf       | brute-force algorithm |
| -kr       | Karp-Rabin algorithm  |
| -lcss     | LCSS algorithm        |

dna_sequence_file: contains a DNA sequence.
pattern_file: contains the pattern that will be searched in the DNA sequence.
second_dna_sequence_file: contains a second DNA sequence.
All files can only contain the following characters:
- c (cytosine)
- g (guanine)
- a (adenine)
- t (thymine)
- \n (new line)

## Examples

```bash
./dnaPatternMatching -bf data/dna1.txt data/pattern1.txt

The pattern was found: 7 times
```

```bash
./dnaPatternMatching -kr data/dna2.txt data/pattern2.txt

The pattern was found: 1 times
```
```bash
./dnaPatternMatching -lcss data/dna1.txt data/dna2.txt

The length of the largest common subsequence is: 241

The distance between the two DNA sequences is: 0.20
```

## Author

Giorgos Argyrides (g.aryrides@outlook.com)