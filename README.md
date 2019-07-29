# LASLO v1.0 - README

# Introduction
Multi-threaded application that searches for hairpin structures in cDNA / mRNA sequences (FASTA or GenBank) with specific consensus sequences in the loop.
It implements RNAfold to predict the structure of the sequence, and UShuffle to generate random sequences with conservation of the k-nucleotide frequency.

# Requirements
JAVA Runtime 8 - Latest update. https://www.java.com/en/download/
It's highly recommended to have a multicore processor to process large sequences.
This version only runs in a Windows x64 SO, because of the external executables. 

# Instructions 

# Credits and Acknowledgments

* RNAfold 2.4.9 - Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler, Ronny Lorenz. https://www.tbi.univie.ac.at/RNA/RNAfold.1.html

* UShuffle - Minghui Jiang, James Anderson, Joel Gillespie, and Martin Mayne. uShuffle: a useful tool for shuffling biological sequences while preserving the k-let counts. BMC Bioinformatics, 9:#192, 2008. 
https://github.com/guma44/ushuffle

* BioJava (Core 4.2.0) - Andreas Prlic; Andrew Yates; Spencer E. Bliven; Peter W. Rose; Julius Jacobsen; Peter V. Troshin; Mark Chapman; Jianjiong Gao; Chuan Hock Koh; Sylvain Foisy; Richard Holland; Gediminas Rimsa; Michael L. Heuer; H. Brandstatter-Muller; Philip E. Bourne; Scooter Willis. BioJava: an open-source framework for bioinformatics in 2012. Bioinformatics (2012) 28 (20): 2693-2695. https://github.com/biojava/biojava

* OpenCSV 4.1 - A Simple CSV Parser for Java under a commercial-friendly Apache 2.0 license  
https://sourceforge.net/projects/opencsv/

# License
 Copyright (C) 2018  David A. Mancilla                                   
 GNU General Public License <http://www.gnu.org/licenses/>. 
