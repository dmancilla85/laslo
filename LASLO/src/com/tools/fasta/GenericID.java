package com.tools.fasta;

/**
 * @author David
 *
 */
public class GenericID extends FastaID {

    final static String HEADER = "SequenceID" + ROW_DELIMITER;

    public GenericID() {
        this.geneID = "";
    }

    @Override
    public String toString() {
        return "GenericID [sequenceID=" + geneID + "]";
    }

    public static String getHeader() {
        return GenericID.HEADER;
    }

    public void setGenericTags(String id) {
        // TODO Auto-generated method stub
        geneID = id.replace(';', '-');
    }

    @Override
    public String toRowCSV() {
        return geneID + ROW_DELIMITER;
    }
}
