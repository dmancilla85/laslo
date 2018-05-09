package com.tools.fasta;

/**
 * @author David
 *
 */
public class FlyBaseFastaID extends FastaID {

    final static String TYPE = "type="; //$NON-NLS-1$
    final static String LOCATION = "loc="; //$NON-NLS-1$
    final static String NAME = "name="; //$NON-NLS-1$
    final static String DBXREF = "dbxref="; //$NON-NLS-1$
    final static String SCORE = "score="; //$NON-NLS-1$
    final static String CHECKSUM = "MD5="; //$NON-NLS-1$
    final static String GENE_ID = "parent="; //$NON-NLS-1$
    final static String RELEASE = "release="; //$NON-NLS-1$
    final static String SPECIES = "species="; //$NON-NLS-1$
    final static String HEADER = "GeneID" + ROW_DELIMITER //$NON-NLS-1$
            + "TranscriptID" + ROW_DELIMITER //$NON-NLS-1$
            + "Name" + ROW_DELIMITER //$NON-NLS-1$
            + "Release" + ROW_DELIMITER //$NON-NLS-1$
            + "DBXREF" + ROW_DELIMITER //$NON-NLS-1$
            + "Score" + ROW_DELIMITER
            + "Checksum" + ROW_DELIMITER //$NON-NLS-1$
            + "Species" + ROW_DELIMITER //$NON-NLS-1$
            ;

    String location;
    String name;
    String dbxref;
    String score;
    String checksum;
    String release;
    String species;

    public FlyBaseFastaID() {
        this.transcriptID = "";
        this.geneID = "";
        this.location = "";     //$NON-NLS-1$
        this.checksum = "";     //$NON-NLS-1$
        this.release = "";      //$NON-NLS-1$
        this.species = "";      //$NON-NLS-1$
        this.name = "";         //$NON-NLS-1$
        this.dbxref = "";       //$NON-NLS-1$
        this.score = "";        //$NON-NLS-1$
    }

    public static String getHeader() {
        return FlyBaseFastaID.HEADER;
    }

    @Override
    public String toString() {
        return "FlyBaseFastaID [transcriptID=" + transcriptID + ", "
                + "geneID=" + geneID + ", location=" + location
                + ", name=" + name + ", dbxref=" + dbxref + ", "
                + "score=" + score + ", checksum=" + checksum
                + ", release=" + release + ", species=" + species + "]";
    }

    public void setFlyBaseTags(String idSequence) {

        int index, index2;
        String idsequence;

        idsequence = idSequence;

        if (idSequence == null || idSequence.length() <= 0) {
            return;
        }

        String aux = idSequence.substring(0, idSequence.indexOf(' '));

        // get Splice number and transcriptID
        index = aux.indexOf('.');

        if (index > 0) {

            // this.id_ensbl.spliceNumber = aux.substring(index + 1).trim();
            transcriptID = aux.substring(0, index - 1).trim();
        } else {
            transcriptID = aux.trim();
        }

        index2 = transcriptID.indexOf(';');

        if (index2 > 0) {
            aux = transcriptID.substring(transcriptID.length(), index2);
        }

        // get location
        index = idSequence.indexOf(FlyBaseFastaID.LOCATION, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.LOCATION.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.LOCATION.length());
            }

            location = aux.trim();
        }

        // get name
        index = idsequence.indexOf(FlyBaseFastaID.NAME, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.NAME.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.NAME.length());
            }

            name = aux.trim();
        }

        // get dbxref
        index = idsequence.indexOf(FlyBaseFastaID.DBXREF, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.DBXREF.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.DBXREF.length());
            }

            dbxref = aux.trim();
        }

        // get score
        index = idsequence.indexOf(FlyBaseFastaID.SCORE, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.SCORE.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.SCORE.length());
            }

            score = aux.trim();
        }

        // get checksum
        index = idsequence.indexOf(FlyBaseFastaID.CHECKSUM, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.CHECKSUM.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.CHECKSUM.length());
            }

            checksum = aux.trim();
        }

        // get GeneID
        index = idsequence.indexOf(FlyBaseFastaID.GENE_ID, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.GENE_ID.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.GENE_ID.length());
            }

            geneID = aux.trim();
        }

        // get release
        index = idsequence.indexOf(FlyBaseFastaID.RELEASE, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.RELEASE.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.RELEASE.length());
            }

            release = aux.trim();
        }

        // get specie
        index = idsequence.indexOf(FlyBaseFastaID.SPECIES, 0);

        if (index > 0) {
            idsequence = idsequence.substring(index);
            index2 = idsequence.indexOf(';');

            if (index2 > 0) {
                aux = idsequence.substring(FlyBaseFastaID.SPECIES.length(), index2);
            } else {
                aux = idsequence.substring(FlyBaseFastaID.SPECIES.length());
            }

            species = aux.trim();
        }
    }

    @Override
    public String toRowCSV() {
        return geneID + ROW_DELIMITER
                + transcriptID + ROW_DELIMITER
                + name + ROW_DELIMITER
                + release + ROW_DELIMITER
                + dbxref + ROW_DELIMITER
                + score + ROW_DELIMITER
                + checksum + ROW_DELIMITER
                + species + ROW_DELIMITER;
    }
}
