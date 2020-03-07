/** 
 *  This program matches patterns between two dna sequences using three different algorithms (brute-force, 
 *  Karp-Rabin, LCSS). Given a dna sequence, a pattern sequence and the name of the algorithm from the user, this 
 *  program executes the selected algorithm to produce information about the pattern matching between the two sequences.
 * 
 *  @author Giorgos Argyrides
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

/* Function prototypes */
int countFileCharacters(char *);
int readFile(char *, char *);
int bruteForce(char *, char *);
int karpRabin(char *, char *);
int lcss(char *, char *);
int hash(char *, int);
int rehash(int, int, int, int);

/** 
 * counts the characters in the file excluding the EOF and the new line characters.
 **/
int countFileCharacters(char *filePath){
	int count = 0;
	char ch;
	FILE *file = NULL;
	if ((file = fopen(filePath, "r")) == NULL){
		printf("\nError: Unable to open %s\n", filePath);
		exit(EXIT_FAILURE);
	}
	while ((ch = fgetc(file)) != EOF){
		if (ch != 'a' && ch != 'c' && ch != 'g' && ch != 't' && ch != '\n'){
			printf("\nerror: a dna sequence can only contain the characters a,c,g,t\n");
			return EXIT_FAILURE;
		}
		if (ch != '\n')
			count++;
	}
	fclose(file);
	return count;
}

/**
 *  Reads all the characters that a dna sequence includes (a,c,g,t) from a file. It returns EXIT_SUCCESS if the 
 *  function completes without errors.
 **/
int readFile(char *sequence, char *filePath){
	FILE *file = NULL;
	if ((file = fopen(filePath, "r")) == NULL){
		printf("\nError: Unable to open %s\n", filePath);
		return EXIT_FAILURE;
	}
	int i = 0;
	char ch = fgetc(file);
	while (ch != EOF){
		// The file may contain \n characters. They will be ignored.
		if (ch != '\n'){
			sequence[i] = ch;
			i++;
		}
		ch = fgetc(file);
	}
	sequence[i] = '\0';
	fclose(file);
	return EXIT_SUCCESS;
}

/** 
 *	Compares the two sequences character by character and prints how many times the pattern sequence is found in the dna
 *  sequence.
 **/
int bruteForce(char *dnaSequence, char *patternSequence){
	int occurrences = 0, i, j, dnaLength = strlen(dnaSequence), patternLength = strlen(patternSequence);
	for (i = 0; i <= dnaLength - patternLength; i++){
		j = 0;
		while (j < patternLength && patternSequence[j] == dnaSequence[i + j]){
			j++;
		}
		if (j == patternLength){
			occurrences++;
		}
	}
	return occurrences;
}

/**  
 *  Implements the Karp-Rabin algorithm for finding subsequences. To avoid multiple character comparison the search is 
 *  divided in two stages: the pre-processing and the actual processing. In the pre-processing stage a hash function is 
 *  used to compute the integer hash value of a dna subsequence (same length with the pattern) and the pattern sequence.
 *  In the processing stage the hash value of the dna subsequence is compared with the hash value of the pattern. The 
 *  next hash value is computed using a rehash function. Returns how many times the pattern sequence is found in the dna 
 *  sequence.
 **/
int karpRabin(char *dnaSequence, char *patternSequence){
	int i, occurrences = 0, dnaLength = strlen(dnaSequence), patternLength = strlen(patternSequence);
	int patHash = hash(patternSequence, patternLength - 1);
	int dnaHash = hash(dnaSequence, patternLength - 1);
	for (i = 0; i <= dnaLength - patternLength; i++){
		if (patHash == dnaHash){
			int j = 0;
			while (j < patternLength && patternSequence[j] == dnaSequence[i + j]){
				j++;
			}
			if (j == patternLength){
				occurrences++;
			}
		}
		dnaHash = rehash((int) dnaSequence[i], dnaHash, (int) dnaSequence[i + patternLength], patternLength - 1);
	}
	return occurrences;
}

/** 
 * Computes the hash value of a character sequence.
 **/
int hash(char *sequence, int length){
	int value = 0, i, p = length;
	for (i = 0; i <= length; i++){
		value += ((int)sequence[i]) * (int)pow(2, p);
		p--;
	}
	value = value % INT_MAX;
	return value;
}

/** 
 * Computes the rehash value of a character sequence.
 **/
int rehash(int before, int hash, int after, int power){
	int value = ((hash - before * (int) pow(2, power)) * 2 + after) % INT_MAX;
	return value;
}

/** 
 * A recursive function for finding the length of the longest common subsequence between two character sequences.
 **/
int lcss(char *sequenceA, char *sequenceB){
	int lcssLength = 0, i, j;
	int lenA = strlen(sequenceA) + 1;
	int lenB = strlen(sequenceB) + 1;
	int *arr = (int *) malloc(lenA * lenB * sizeof(int)); 
	for(i = 0; i < lenA; i++){
		*(arr + i * lenB) = 0;
	}
	for(j = 0; j < lenB; j++){
		*(arr + j)  = 0;
	}
	for(i = 1; i < lenA; i++){
		for(j = 1; j < lenB; j++){
			if(sequenceA[i-1] == sequenceB[j-1]){
				*(arr + i * lenB + j) = *(arr + (i - 1) * lenB + j - 1) + 1;
			}
			else{
				*(arr + i * lenB + j) = max(*(arr + i * lenB + j - 1), *(arr + (i - 1) * lenB + j));
			}
		}
	}
	lcssLength = *(arr + (lenA - 1) * lenB + lenB - 1);
	free(arr);
	return lcssLength;
}

int main(int argc, char **argv){	
	if (argc != 4){
		printf("\nError: invalid arguments.\nUsage:\n ./dnaPatternMatching <algorithm> <dna_sequence_file> "
		"<pattern_file or second_dna_sequence_file>\n");
		exit(EXIT_FAILURE);
	}
	char *dnaSequence = NULL;
	int dnaLength = countFileCharacters(argv[2]);
	dnaSequence = (char *) malloc(dnaLength * sizeof(char));
	if (dnaSequence == NULL){
		printf("\nError: failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}
	char *patternSequence = NULL;
	int patternLength = countFileCharacters(argv[3]);
	patternSequence = (char *) malloc(patternLength * sizeof(char));
	if (dnaSequence == NULL){
		printf("\nError: failed to allocate memory.\n");
		exit(EXIT_FAILURE);
	}
	if (readFile(dnaSequence, argv[2]) == EXIT_FAILURE){
		exit(EXIT_FAILURE);
	}
	if (readFile(patternSequence, argv[3]) == EXIT_FAILURE){
		exit(EXIT_FAILURE);
	}
	if (dnaLength < patternLength && (strcmp(argv[1], "-lcss") != 0)){
		printf("\nError: the dna sequence must contain more characters than the pattern sequence\n");
		exit(EXIT_FAILURE);
	}
	if (strcmp(argv[1], "-bf") == 0){
		int occurrences = bruteForce(dnaSequence, patternSequence);
		printf("\nThe pattern was found: %d times\n\n", occurrences);
	}
	else if (strcmp(argv[1], "-kr") == 0){
		int occurrences = karpRabin(dnaSequence, patternSequence);
		printf("\nThe pattern was found: %d times\n\n", occurrences);
	}
	else if (strcmp(argv[1], "-lcss") == 0){
		int lcssLength = lcss(dnaSequence, patternSequence);
		printf("\nThe length of the largest common subsequence is: %d\n", lcssLength);
		float distance = 1 - (float)lcssLength / min(strlen(dnaSequence), strlen(patternSequence));
		printf("\nThe distance between the two DNA sequences is: %.2f\n", distance);
	}
	else{
		printf("\nError: invalid arguments.\nUsage:\n ./dnaPatternMatching <algorithm> <dna_sequence_file> "
		"<pattern_file or second_dna_sequence_file>\n");
		exit(EXIT_FAILURE);
	}
	// free(dnaSequence);
	// free(patternSequence);
	return 0;
}
