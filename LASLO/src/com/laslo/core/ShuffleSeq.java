/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.laslo.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import static java.lang.System.out;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author David A. Mancilla
 */
public class ShuffleSeq {

    public static String commandShuffleSeq = "./shuffleseq.exe";
    protected static String randDir = "/shuffled";

    public static void main(String[] args) {
  
            Runtime rt = Runtime.getRuntime();
            String path = "C:/Users/David/Desktop/profile/";
            String ext = ".fasta";
            String fileName = "repressed_non_bound";
            
            makeShuffleSequences(path, fileName, ext, 1000);

    }

    /**
     * 
     * @param path
     * @param fileName
     * @param extension
     * @param nRandoms 
     */
    public static void makeShuffleSequences(String path, String fileName,
            String extension, int nRandoms) {
        Runtime rt;
        rt = Runtime.getRuntime();
        String origin = path + fileName + extension;
        String destiny = "";
        boolean mkdirs = new File(path + randDir).mkdirs();
        
        for (int i = 1; i <= nRandoms; i++) {
            try {

                destiny = path + randDir + '/' + fileName + "_rnd" + i + extension;

                if (new File(destiny).exists()) {
                    boolean delete = (new File(destiny)).delete();
                }

                Process pr = rt.exec(commandShuffleSeq + " -sequence "
                        + origin + " -outseq " + destiny);
                int exitVal = pr.waitFor();
                System.out.println("From " + origin + " to " + destiny);
                System.out.println("System exit code: " + exitVal);
            } catch (IOException | InterruptedException ex) {
                Logger.getLogger(ShuffleSeq.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }
}
