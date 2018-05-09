/**
 *
 */
package com.tools.fasta;

/**
 * @author David
 *
 */
public class EnsemblFastaID extends FastaID {

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

    public EnsemblFastaID() {
        this.transcriptID = ""; //$NON-NLS-1$
        this.geneID = ""; //$NON-NLS-1$
        this.geneBiotype = ""; //$NON-NLS-1$
        this.transcriptBiotype = ""; //$NON-NLS-1$
        this.geneSymbol = ""; //$NON-NLS-1$
        this.description = ""; //$NON-NLS-1$
        this.spliceNumber = ""; //$NON-NLS-1$

    }

    public String getGeneBiotype() {
        return geneBiotype;
    }

    public void setGeneBiotype(String geneBiotype) {
        this.geneBiotype = geneBiotype;
    }

    public String getTranscriptBiotype() {
        return transcriptBiotype;
    }

    public void setTranscriptBiotype(String transcriptBiotype) {
        this.transcriptBiotype = transcriptBiotype;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getSpliceNumber() {
        return spliceNumber;
    }

    public void setSpliceNumber(String spliceNumber) {
        this.spliceNumber = spliceNumber;
    }

    @Override
    public String toString() {
        return "EnsemblFastaID [transcriptID=" + transcriptID + ", geneID=" + geneID + ", geneBiotype=" + geneBiotype //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
                + ", trancriptBiotype=" + transcriptBiotype + ", geneSymbol=" + geneSymbol + ", description=" //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
                + description + ", spliceNumber=" + spliceNumber + "]"; //$NON-NLS-1$ //$NON-NLS-2$
    }

    public static String getHeader() {
        return EnsemblFastaID.HEADER;
    }

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
