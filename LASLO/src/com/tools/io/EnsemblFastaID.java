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
public class EnsemblFastaID extends SourceFile {

    public final String GENE = "gene:"; //$NON-NLS-1$
    public final static String GENE_BIOTYPE = "gene_biotype:"; //$NON-NLS-1$
    public final static String TRANSCRIPT_BIOTYPE = "transcript_biotype:"; //$NON-NLS-1$
    public final static String GENE_SYMBOL = "gene_symbol:"; //$NON-NLS-1$
    public final static String DESCRIPTION = "description:"; //$NON-NLS-1$

    final static String HEADER = "GeneID" + ROW_DELIMITER //$NON-NLS-1$
            + "GeneSymbol" + ROW_DELIMITER //$NON-NLS-1$
            + "GeneBioType" + ROW_DELIMITER //$NON-NLS-1$
            + "TranscriptID" + ROW_DELIMITER //$NON-NLS-1$
            + "TranscriptBiotype" + ROW_DELIMITER //$NON-NLS-1$
            + "#Splice" + ROW_DELIMITER //$NON-NLS-1$
            + "Description" + ROW_DELIMITER //$NON-NLS-1$
            ;

    protected String geneBiotype;
    protected String transcriptBiotype;
    protected String description;
    protected String spliceNumber;

    /**
     * 
     */
    public EnsemblFastaID() {
        this.transcriptID = ""; //$NON-NLS-1$
        this.geneID = ""; //$NON-NLS-1$
        this.geneBiotype = ""; //$NON-NLS-1$
        this.transcriptBiotype = ""; //$NON-NLS-1$
        this.geneSymbol = ""; //$NON-NLS-1$
        this.description = ""; //$NON-NLS-1$
        this.spliceNumber = ""; //$NON-NLS-1$

    }

    /**
     * 
     * @return 
     */
    public String getGeneBiotype() {
        return geneBiotype;
    }

    /**
     * 
     * @param geneBiotype 
     */
    public void setGeneBiotype(String geneBiotype) {
        this.geneBiotype = geneBiotype;
    }

    /**
     * 
     * @return 
     */
    public String getTranscriptBiotype() {
        return transcriptBiotype;
    }

    /**
     * 
     * @param transcriptBiotype 
     */
    public void setTranscriptBiotype(String transcriptBiotype) {
        this.transcriptBiotype = transcriptBiotype;
    }

    /**
     * 
     * @return 
     */
    public String getDescription() {
        return description;
    }

    /**
     * 
     * @param description 
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * 
     * @return 
     */
    public String getSpliceNumber() {
        return spliceNumber;
    }

    /**
     * 
     * @param spliceNumber 
     */
    public void setSpliceNumber(String spliceNumber) {
        this.spliceNumber = spliceNumber;
    }

    /**
     * 
     * @return 
     */
    @Override
    public String toString() {
        return "EnsemblFastaID [transcriptID=" + transcriptID + ", geneID=" + geneID + ", geneBiotype=" + geneBiotype //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
                + ", trancriptBiotype=" + transcriptBiotype + ", geneSymbol=" + geneSymbol + ", description=" //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
                + description + ", spliceNumber=" + spliceNumber + "]"; //$NON-NLS-1$ //$NON-NLS-2$
    }

    /**
     * 
     * @return 
     */
    public static String getHeader() {
        return EnsemblFastaID.HEADER;
    }

    /**
     * 
     * @param idSequence 
     */
    public void setEnsemblTags(String idSequence) {

        int index, index2;
        String idsequence;
        String aux;

        if (idSequence == null || idSequence.length() <= 0) {
            return;
        }

        if (idSequence.indexOf(' ') < 0) {
            setGeneID(idSequence);
            return;
        }

        aux = idSequence.substring(0, idSequence.indexOf(' '));

        // Just to don't break the .csv file
        idsequence = idSequence.replaceAll(";", ":"); //$NON-NLS-1$ //$NON-NLS-2$

        // get Splice number and transcriptID
        index = aux.indexOf('.');

        if (index > 0) {
            setSpliceNumber(aux.substring(index + 1).trim());
            setTranscriptID(aux.substring(0, index - 1).trim());
        } else {
            setTranscriptID(aux.trim());
        }

        // get GeneID
        index = idSequence.indexOf(GENE, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(' ');

            if (index2 > 0) {
                aux = idsequence.substring(GENE.length(), index2);
            } else {
                aux = idsequence.substring(GENE.length());
            }

            setGeneID(aux.trim());
        }

        // get Gene Biotype
        index = idsequence.indexOf(GENE_BIOTYPE, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(' ');

            if (index2 > 0) {
                aux = idsequence.substring(GENE_BIOTYPE.length(), index2);
            } else {
                aux = idsequence.substring(GENE_BIOTYPE.length());
            }

            setGeneBiotype(aux.trim());
        }
        // get Trancript Biotype
        index = idsequence.indexOf(EnsemblFastaID.TRANSCRIPT_BIOTYPE, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(' ');

            if (index2 > 0) {
                aux = idsequence.substring(EnsemblFastaID.TRANSCRIPT_BIOTYPE.length(), index2);
            } else {
                aux = idsequence.substring(EnsemblFastaID.TRANSCRIPT_BIOTYPE.length());
            }

            setTranscriptBiotype(aux.trim());
        }

        // get Gene Symbol
        index = idsequence.indexOf(EnsemblFastaID.GENE_SYMBOL, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(' ');

            if (index2 > 0) {
                aux = idsequence.substring(EnsemblFastaID.GENE_SYMBOL.length(), index2);
            } else {
                aux = idsequence.substring(EnsemblFastaID.GENE_SYMBOL.length());
            }

            setGeneSymbol(aux.trim());
        }

        // Description
        index = idsequence.indexOf(EnsemblFastaID.DESCRIPTION, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.lastIndexOf(':');

            if (index2 > 0) {
                aux = idsequence.substring(EnsemblFastaID.DESCRIPTION.length(), idsequence.lastIndexOf(' '));
            } else {
                aux = idsequence.substring(EnsemblFastaID.DESCRIPTION.length());
            }

            setDescription(aux.trim());
        }

    }

    /**
     * 
     * @return 
     */
    @Override
    public String toRowCSV() {
        return geneID + ROW_DELIMITER
                + geneSymbol + ROW_DELIMITER
                + geneBiotype + ROW_DELIMITER
                + transcriptID + ROW_DELIMITER
                + transcriptBiotype + ROW_DELIMITER
                + spliceNumber + ROW_DELIMITER
                + description + ROW_DELIMITER;
    }
}
