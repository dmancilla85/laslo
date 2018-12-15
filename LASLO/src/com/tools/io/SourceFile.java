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
package com.tools.io;

/**
 * @author David
 *
 */
public class SourceFile {

    final static String FIELD_DELIMITER = "|"; //$NON-NLS-1$

    /**
     *
     */
    public final static String ROW_DELIMITER = ";"; //$NON-NLS-1$

    /**
     *
     */
    protected String transcriptID;

    /**
     *
     */
    protected String geneID;

    /**
     *
     */
    protected String geneSymbol;

    /**
     *
     */
    protected String fieldDelimiter;

    /**
     *
     */
    protected String rowDelimiter;

    /**
     *
     */
    public static final String LOG_EXT = ".log";

    /**
     *
     */
    public static final String CSV_EXT = ".csv";

    /**
     *
     */
    public static final String FASTA_EXT = ".fasta";

    /**
     *
     */
    public static final String FASTA_EXT_2 = ".fa";

    /**
     *
     */
    public static final String GENBANK_EXT = ".gb";

    /**
     *
     * @param header
     * @return
     */
    public static InputSequence detectHeader(String header) {

        InputSequence out = InputSequence.GENERIC;

        if(header.contains("genbank")){
            return InputSequence.GENBANK;
        }
        
        if(header.contains("@")){
            return InputSequence.GENBANK;
        }
        
        if (header.contains("gene:")
                && header.contains("gene_biotype:")
                && header.contains("transcript_biotype:")
                && header.contains("gene_symbol:")
                && header.contains("description:")) {
            out = InputSequence.ENSEMBL;
        } else {
            if (header.contains("type=")
                    && header.contains("loc=")
                    && header.contains("name=")
                    && header.contains("dbxref=")
                    && header.contains("score=")
                    && header.contains("MD5=")
                    && header.contains("release=")) {
                out = InputSequence.FLYBASE;
            } else {
                if (header.contains("|") || header.contains("=")) {
                    out = InputSequence.BIOMART;
                }
            }
        }

        return out;
    }

    /**
     *
     * @return
     */
    public String getTranscriptID() {
        return transcriptID;
    }

    /**
     *
     * @param transcriptID
     */
    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }

    /**
     *
     * @return
     */
    public String getGeneID() {
        return geneID;
    }

    /**
     *
     * @param geneID
     */
    public void setGeneID(String geneID) {
        this.geneID = geneID;
    }

    /**
     *
     * @return
     */
    public String getGeneSymbol() {
        return geneSymbol;
    }

    /**
     *
     * @param geneSymbol
     */
    public void setGeneSymbol(String geneSymbol) {
        this.geneSymbol = geneSymbol;
    }

    /**
     *
     * @return
     */
    public String getRowDelimiter() {
        return rowDelimiter;
    }

    /**
     *
     * @param rowDelimiter
     */
    public void setRowDelimiter(String rowDelimiter) {
        this.rowDelimiter = rowDelimiter;
    }

    /**
     *
     * @return
     */
    public String getFieldDelimiter() {
        return fieldDelimiter;
    }

    /**
     *
     * @param fieldDelimiter
     */
    public void setFieldDelimiter(String fieldDelimiter) {
        this.fieldDelimiter = fieldDelimiter;
    }

    /**
     *
     */
    public SourceFile() {
        this.fieldDelimiter = FIELD_DELIMITER;
        this.rowDelimiter = ROW_DELIMITER;
    }

    /**
     *
     * @return
     */
    public static String getHeader() {
        return "NO HEADER";
    }

    /**
     *
     * @return
     */
    @Override
    public String toString() {
        return "geneID = " + geneID + " - "
                + "transcriptID = " + transcriptID;
    }

    /**
     *
     * @return
     */
    public String toRowCSV() {
        return geneID + ROW_DELIMITER
                + transcriptID + ROW_DELIMITER;
    }
}
