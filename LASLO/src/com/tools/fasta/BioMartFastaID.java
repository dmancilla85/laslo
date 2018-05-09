/**
 *
 */
package com.tools.fasta;

/**
 * @author David A. Mancilla
 *
 */
public class BioMartFastaID extends FastaID {

    protected static int nCols = 6;

    /**
     *
     */
    protected static String columns[] = new String[nCols];

    public static int getNumColumns() {
        return columns.length;
    }

    public static void setColumns(String cols[]) {
        columns = cols;
    }

    public static String getHeader() {
        String header = "";
        
        for (int i = 0; i < nCols; i++) {
            header += "Column" + (i + 1) + ROW_DELIMITER;
        }

        return header;
    }

    public void setBioMartTags(String idSequence) {

        String auxSeq = idSequence.replace('=', '@');
        auxSeq = auxSeq.replace(';', '-');
        String [] cols = auxSeq.split("@");
        
        for(int i = 0; i < cols.length && i < 6; i++){
            columns[i] = cols[i];
        }
    }

    @Override
    public String toRowCSV() {
        String row = "";
        int i;

        for (i = 0; i < columns.length; i++) {
            row += columns[i] + ROW_DELIMITER;
        }

        return row;
    }
}
