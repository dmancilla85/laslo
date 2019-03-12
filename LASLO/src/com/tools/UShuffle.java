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
package com.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.features.TextFeature;
import static org.biojava.nbio.core.sequence.io.FastaReaderHelper.readFastaDNASequence;

/**
 *
 * @author David A. Mancilla
 */
public class UShuffle {

    private static final String COMMAND_SHUFFLE = "./ext/ushuffle.exe -s ";
    private static final String RANDOM_PATH = "\\shuffled";

    /**
     *
     * @param args
     */
    public static void main(String[] args) {

        String path = "C:\\Users\\David\\Documents\\NetBeans"
                + "Projects\\lasloProject\\LASLO\\ext";
        String ext = ".fasta";
        String fileName = "fruitfly_chen";

        System.out.println("Iniciando...");
        // Generation of the iterator of {id,sequence}
        LinkedHashMap<String, DNASequence> fasta = null;
        try {
            fasta = readFastaDNASequence(new File("./ext/" + fileName + ext), false);
        } catch (IOException ex) {
            Logger.getLogger(UShuffle.class.getName()).log(Level.SEVERE, null, ex);
        }
        makeShuffleSequences(path, fileName, fasta, 5, 2, false);

    }

    /**
     *
     * @return
     */
    public static String getRandomDir() {
        return RANDOM_PATH;
    }

    /**
     *
     * @param path      Path of the file to shuffle
     * @param filename  Name of the file to shuffle
     * @param fasta     Hashmap of DNASequence's from BioJava
     * @param nRandoms  Number of random sequences to generate
     * @param k         Value of k-let permutations
     * @param isGenBank It tells if it's a GenBank file
     * @return          Last process exit value (an integer)
     */
    public static int makeShuffleSequences(String path, String filename,
            LinkedHashMap<String, DNASequence> fasta, int nRandoms, int k,
            boolean isGenBank) {
        Runtime rt;
        char sep = '@';
        int exitVal = 0;
        String aux = "", gene, synonym, note;
        rt = Runtime.getRuntime();
        String fileNameWithOutExt = filename.replaceFirst("[.][^.]+$", "");
        String destiny = "", sequence = "", header = "", id = "", cds = "";
        boolean mkdirs, delete;
        mkdirs = new File(path + RANDOM_PATH).mkdirs();

        for (int i = 1; i <= nRandoms; i++) {
            try {

                destiny = path + RANDOM_PATH + '\\' + fileNameWithOutExt
                        + "_rnd" + i + ".fa";
                File file = new File(destiny);

                if (file.exists()) {
                    delete = (new File(destiny)).delete();
                }

                file.createNewFile();

                FileWriter fw = new FileWriter(file.getAbsoluteFile(), true);
                BufferedWriter bw = new BufferedWriter(fw);
                int j = 1;
                for (Map.Entry<String, DNASequence> entry : fasta.entrySet()) {
                    DNASequence element = entry.getValue();

                    if (!isGenBank) {
                        header = element.getOriginalHeader();
                    } else {
                        Map qual = ((TextFeature) element.getFeaturesByType("gene")
                                .toArray()[0]).getQualifiers();

                        if (!qual.isEmpty()) {
                            gene = ((Qualifier) ((ArrayList) (qual.get("gene"))).get(0))
                                    .getValue();
                            synonym = ((Qualifier) ((ArrayList) (qual.get("gene_synonym"))).get(0))
                                    .getValue();
                            note = ((Qualifier) ((ArrayList) (qual.get("note"))).get(0))
                                    .getValue();
                            id = element.getAccession().getID();
                            cds = ((TextFeature) element.getFeaturesByType("CDS")
                                    .toArray()[0]).getSource();

                            header = gene + sep + synonym + sep
                                    + note + sep + id + sep + cds;
                        }
                    }

                    sequence = element.getSequenceAsString();
                    String cmd = COMMAND_SHUFFLE + sequence
                            + " -n 1 -k " + k;
                    Process pr = rt.exec(cmd);

                    try (InputStream in = pr.getInputStream()) {
                        int c;

                        while ((c = in.read()) != -1) {

                            if (c != '\n') {
                                aux += (char) c;
                            } else {
                                bw.write(">" + header);
                                bw.newLine();
                                bw.write(aux);
                                bw.newLine();
                                bw.newLine();
                                aux = "";
                            }

                        }
                    }

                    exitVal = pr.waitFor();

                    /*if(exitVal != 0) {
                    System.out.println(destiny + ": Sequence " + j++
                        + " - ["+ exitVal + "]");
                    }*/
                }

                bw.close();
                bw = null;
                fw = null;

            } catch (IOException | InterruptedException ex) {
                System.out.println("Error: " + ex.getMessage());
            }
        }
        
        return exitVal;
    }
}
