/*
 * Copyright (C) 2018 David A. Mancilla
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package com.laslo.core;

/**
 * @author David
 *
 */
public class PairmentAnalizer {

    /**
     * 
     * @param base1
     * @param base2
     * @return 
     */
    public static boolean isComplementaryDNA(char base1, char base2) {
        // ok
        boolean isComplement = false;

        switch (base1) {
            case 'A':
                isComplement = (base2 == 'T');
                break;
            case 'T':
                isComplement = (base2 == 'A');
                break;
            case 'G':
                isComplement = (base2 == 'C');
                break;
            case 'C':
                isComplement = (base2 == 'G');
                break;
        }

        return isComplement;
    }

    /**
     * 
     * @param source
     * @return 
     */
    public static String reverseIt(String source) {
        int i, len = source.length();
        StringBuilder dest = new StringBuilder(len);

        for (i = (len - 1); i >= 0; i--) {
            dest.append(source.charAt(i));
        }

        return dest.toString();
    }

    /**
     * 
     * @param base1
     * @param base2
     * @return 
     */
    public static boolean isComplementaryRNA(char base1, char base2) {
        // ok
        boolean isComplement = false;

        /**
         * The nucleic acid codes are: A --> adenosine M --> A C (amino) C -->
         * cytidine S --> G C (strong) G --> guanine W --> A T (weak) T -->
         * thymidine B --> G T C U --> uridine D --> G A T R --> G A (purine) H
         * --> A C T Y --> T C (pyrimidine) V --> G C A K --> G T (keto) N --> A
         * G C T (any) - gap of indeterminate length
         */
        switch (base1) {
            case 'A':
                isComplement = (base2 == 'U') || (base2 == 'W')
                        || (base2 == 'K') || (base2 == 'Y')
                        || (base2 == 'B') || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'U':
                isComplement = (base2 == 'A') || (base2 == 'W')
                        || (base2 == 'M') || (base2 == 'R')
                        || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'G':
                isComplement = (base2 == 'C') || (base2 == 'S')
                        || (base2 == 'M') || (base2 == 'Y')
                        || (base2 == 'B')
                        || (base2 == 'H') || (base2 == 'V');
                break;
            case 'C':
                isComplement = (base2 == 'G') || (base2 == 'S')
                        || (base2 == 'K') || (base2 == 'R')
                        || (base2 == 'B');
                break;
        }

        return isComplement;
    }

    /**
     * 
     * @param base1
     * @param base2
     * @return 
     */
    public static boolean isMismatch(char base1, char base2) {
        boolean isMismatch;
        isMismatch = !isComplementaryDNA(base1, base2)
                && !isComplementaryRNA(base1, base2)
                && !isComplementaryRNAWooble(base1, base2);

        return isMismatch;
    }

    /**
     * The nucleic acid codes are: 
     * A --> adenosine 
     * M --> A C (amino) 
     * C --> cytidine 
     * S --> G C (strong) 
     * G --> guanine 
     * W --> A T (weak) 
     * T --> thymidine 
     * B --> G T C 
     * U --> uridine 
     * D --> G A T 
     * R --> G A (purine) 
     * H --> A C T 
     * Y --> T C (pyrimidine) 
     * V --> G C A 
     * K --> G T (keto) 
     * N --> A G C T (any) 
     * - gap of indeterminate length
     * @param base1
     * @param base2
     * @return 
     */
    public static boolean isComplementaryRNAWooble(char base1, char base2) {

        boolean isComplement = false;

        switch (base1) {
            case 'U':
                isComplement = (base2 == 'G') || (base2 == 'S')
                        || (base2 == 'K') || (base2 == 'R')
                        || (base2 == 'B') || (base2 == 'I');
                break;

            case 'G':
                isComplement = (base2 == 'U') || (base2 == 'W')
                        || (base2 == 'K') || (base2 == 'Y')
                        || (base2 == 'B') || (base2 == 'D')
                        || (base2 == 'H') || (base2 == 'V');
                break;

            case 'I':
                isComplement = (base2 == 'U') || (base2 == 'A')
                        || (base2 == 'C');
                break;

            case 'A':
                isComplement = (base2 == 'I');
                break;

            case 'C':
                isComplement = (base2 == 'I');
                break;

        }

        return isComplement;
    }

    /**
     * 
     * @param dnaSequence
     * @return 
     */
    public static String toComplementaryRNA(String dnaSequence) {
        String rnaSequence = ""; //$NON-NLS-1$
        int i;

        for (i = 0; i < dnaSequence.length(); i++) {
            switch (dnaSequence.charAt(i)) {
                case 'A':
                    rnaSequence = rnaSequence.concat("U"); //$NON-NLS-1$
                    break;
                case 'T':
                    rnaSequence = rnaSequence.concat("A"); //$NON-NLS-1$
                    break;
                case 'G':
                    rnaSequence = rnaSequence.concat("C"); //$NON-NLS-1$
                    break;
                case 'C':
                    rnaSequence = rnaSequence.concat("G"); //$NON-NLS-1$
                    break;
                default:
                    rnaSequence = rnaSequence.concat(String.valueOf(
                            dnaSequence.charAt(i))); // test
            }
        }

        return rnaSequence;
    }

    /**
     * 
     * @param fastaPattern
     * @return 
     */
    public static String toRegularExpression(String fastaPattern) {
        String regExp = null;
        boolean extendedReplace = true;

        regExp = fastaPattern.replaceAll("N", "[AUTGC]"); //$NON-NLS-1$ 

        if (extendedReplace) {
            regExp = regExp.replaceAll("W", "[AUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("S", "[CG]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("M", "[AC]"); //$NON-NLS-1$
            regExp = regExp.replaceAll("K", "[GUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("R", "[GA]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("Y", "[CUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("B", "[GCUT]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("D", "[GUTA]"); //$NON-NLS-1$
            regExp = regExp.replaceAll("H", "[CUTA]"); //$NON-NLS-1$ 
            regExp = regExp.replaceAll("V", "[CGA]"); //$NON-NLS-1$ 
        }
        return regExp;
    }

    /**
     * 
     * @param rna
     * @param sequence
     * @param nMin
     * @return 
     */
    public static int findSlippageSequence(String rna, String sequence, 
            int nMin) {
        /**
         * To search for repeated tracts as CA[N]
         */
        int init = 0, n = 0, position = 0, indexOf;

        indexOf = rna.indexOf(sequence, init);

        while (indexOf > 0) {

            if (init == 0) {
                position = indexOf;
            }

            init = indexOf + sequence.length();
            n++;

            indexOf = rna.indexOf(sequence, init);
        }

        if (n >= nMin) {
            return position;
        } else {
            return 0;
        }
    }

    /**
     * 
     * @param loop
     * @return 
     */
    public static boolean checkInternalPairments(String loop) {
        char base1, base2;
        boolean ret = true;

        // Check if loop is good (There's no pairment between his bases)
        if (loop.length() > 5) {
            for (int i = 0; i < 1/*loop.length()/2 - 1*/; i++) {

                base1 = loop.charAt(i);
                base2 = loop.charAt(loop.length() - 1 - i);

                if (PairmentAnalizer.isComplementaryRNAWooble(base1, base2) || 
                        PairmentAnalizer.isComplementaryRNA(base1, base2)) {
                    ret = false;
                }
            }
        }
        return ret;
    }

}
