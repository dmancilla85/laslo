# LASLO v1.01 - README

## Introduction

Multi-threaded application that searches for hairpin structures in cDNA / mRNA sequences (FASTA or GenBank) with specific consensus sequences in the loop. It implements RNAfold to predict the structure of the sequence, and UShuffle to generate random sequences with conservation of the k-nucleotide frequency.

## Requirements

JAVA Runtime 8 - Latest update. https://www.java.com/en/download/ It's highly recommended to have a multicore processor to process large sequences. This version only runs in a Windows x64 SO, because of the external executables.

## Latest release

- 25/8. Fixed problem with generation of multiple shuffled sequences.

## Instructions

> **Input modes**
> The File input tab allows you to select a FASTA file, or a folder containing various FASTA files to process.
> The NCBI input tab allows you to enter a list of NCBI id accession (separated by comma) to download the GenBank files from NCBI (e.g. NM_130489), the ID accession must be of a transcript sequence. The GenBank obtained will be downloaded to the destination folder.

> **Parameters**
> You must select a destination folder to save the results and the shuffled files generated and/or downloaded from NCBI.
>
>> **Shuffled sequences**
>> If you check the option to generate shuffled sequences, you must select the number of shuffled sequences to make, and the k-let parameter. 
>> The parameter k indicates the number of nucleotides that will be taken together when shuffling the sequence. The default value is k=2, this value conserves the stacking energy between bases. 
>
>> **Fold mode**
>> If you select the 'full fold' mode, the entire sequence will be processed with the RNAFold and then, the application will look for stem-loops matching the loop patterns ingresed. If you do not select the full fold mode, only one section of the entire sequence will be folded with each matching loop. This section will contain the pattern sought in the center and at the ends, two frames with the same or less number of bases as the maximum length indicated for the stem.
>
>> **Temperature**
>> Sets the folding temperature to RNAfold. By default 25ºC.
>
>> **Stem length**
>> Determines the range desired of the stem length. 
>
>> **Loop patterns**
>> Here you enter the loop patterns to find (using the IUPAC fasta nomenclature). If you enter a list, must be separated by comma.
>> *IUPAC nucleotide codes: https://www.bioinformatics.org/sms/iupac.html*
>
>> **Search in anti-sense too**
>> Find the pattern ingresed in both senses.
>
>> **Find a pattern anywhere**
>> Find and additional simple pattern in the entire sequence, the output indicate the positions and number of ocurrences.

> **Output table**		
>> ***Note 1. The first columns may vary, it depends on the anotations in the FASTA files.***
>> ***Actually is able to parse Ensembl, FlyBase, GenBank and BIOMart headers. Each of these formats have different initial columns.***

>>	**LoopPattern** (nominal): The pattern matched by the current loop.

>>	**Sense**: (+) sense or (-) nonsense. Useful if the "Search in anti-sense" is checked.

>>	**TerminalPair** (nominal):​ Loop closing pairs. In 5' to 3' order.

>>	**N-2** (char): Precursor nucleotide to the 5 'loop predecessor. (Loop start 5' minus two nucleotide).

>>	**N-1** (char)​: Predecesor nucleotide of the loop in 5'.

>>	**N2** (char)​: Nucleotide in the 2nd position of the loop.

>>	**N5** (char): Nucleotide in the 5th position of the loop.

>>	**N6** (char): Nucleotide in the 6th position of the loop.

>>	**N7** (char)​: Nucleotide in the 7th position of the loop.

>>	**N8** (char)​: Nucleotide in the 8th position of the loop.

>> 	***Note 2: Originally intended for patterns like CNGGNNNN. It will be generalized soon for al N variables.***
>>	**Loop** (string): Loop sequence.

>>	**StemLoopSequence** (string) Stem-loop full sequence.

>>	**Additional5Seq** (string​): A sequence in 5'.

>>	**Additional3Seq​** (string)​: A sequence in 3'.

>>	**PredictedStructure** (string): Vienna bracket structure with asterisks in the wooble sites.

>>	**ViennaBracketStr** (string): Vienna bracket structure returned by RNAFold.

>>	**Pairments** (integer): Count of pairments in the stem.

>>	**WooblePairs** (integer): Count of wooble pairs in the stem.

>>	**Bulges** (integer​): Count of bulges in the stem-loop.

>>	**InternalLoops** (integer​): Count of internal loops in the stem-loop.

>>	**SequenceLength** (integer​): Length (bases) of the full sequence.

>>	**StartsAt** (numeric​): Start position (first base) of the stem-loop (base number starting from the beginning).

>>	**EndsAt** (numeric​): End position of the stem-loop (last base).

>>	**A_PercentSequence** (double​): Percent of A in the sequence analized.

>>	**C_PercentSequence** (double​): Percent of C in the sequence analized.

>>	**G_PercentSequence** ​(double​): Percent of G in the sequence analized.

>>	**U_PercentSequence**​ (double​): Percent of U in the sequence analized.

>>	**AU_PercentPairs**​	(double​): Percent of AU pairments in the stem-loop.

>>	**CG_PercentPairs**​	(double​): Percent of CG pairments in the stem-loop.

>>	**GU_PercentPairs**​	(double​): Percent of GU pairments in the stem-loop.

>>	**PurinePercentStem** (double​): Percent of A/G bases inside the stem-loop.

>>	**RnaFoldMFE** (double): Minimun free energy obtained from RNAFold for the sequence (kcal/mol).

>>	**RelativePosition** (double): Start position of the stem-loop measured as StartsAt/SequenceLength.

>>	**AdditionalSeqMatches** (integer)​: Count of matches for the additional sequence (If it applies).

>>	**AdditionalSeqPosition**​ (list): List of start positions of the additional sequence coincidences.

## Credits and Acknowledgments

1. RNAfold 2.4.9 - Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler, Ronny Lorenz. https://www.tbi.univie.ac.at/RNA/RNAfold.1.html

2. UShuffle - Minghui Jiang, James Anderson, Joel Gillespie, and Martin Mayne. uShuffle: a useful tool for shuffling biological sequences while preserving the k-let counts. BMC Bioinformatics, 9:#192, 2008. https://github.com/guma44/ushuffle

3. BioJava (Core 4.2.0) - Andreas Prlic; Andrew Yates; Spencer E. Bliven; Peter W. Rose; Julius Jacobsen; Peter V. Troshin; Mark Chapman; Jianjiong Gao; Chuan Hock Koh; Sylvain Foisy; Richard Holland; Gediminas Rimsa; Michael L. Heuer; H. Brandstatter-Muller; Philip E. Bourne; Scooter Willis. BioJava: an open-source framework for bioinformatics in 2012. Bioinformatics (2012) 28 (20): 2693-2695. https://github.com/biojava/biojava

4. OpenCSV 4.1 - A Simple CSV Parser for Java under a commercial-friendly Apache 2.0 license https://sourceforge.net/projects/opencsv/

5. Icon [Folder](https://thenounproject.com/icon/53223/) by [Melissa Holterman] (https://thenounproject.com/swiffermuis/) is licensed under [CC BY 4.0] (https://creativecommons.org/licenses/by/4.0/)

6. Icon [DNA](https://thenounproject.com/search/?q=dna&i=1088243) by [Creative Stall] (https://thenounproject.com/creativestall/) is licensed under [CC BY 4.0] (https://creativecommons.org/licenses/by/4.0/)

## License

*Copyright (C) 2018 David A. Mancilla GNU General Public License http://www.gnu.org/licenses/.*
