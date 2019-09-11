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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import static java.lang.System.out;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author David A. Mancilla
 */
public class RNAfold {

    private final static String COMMAND_RNAFOLD = "./ext/RNAfold.exe";
    //private final static String RNAFOLD_ARGS = " -d2 --noLP --noPS";
    private String structure;
    private double mfe;
    private Double temperature;
    private boolean avoidLonelyPairs;
    private Exception exc;

    /**
     *
     */
    public RNAfold() {
        this.mfe = 0.0;
        this.structure = "";
        this.temperature = 25.00;
        this.avoidLonelyPairs = true;
        this.exc = null;
    }

    /**
     *
     * @param temp
     * @param avoidLonelyPairs
     */
    public RNAfold(double temp, boolean avoidLonelyPairs) {
        this.mfe = 0.0;
        this.structure = "";
        this.temperature = temp;
        this.avoidLonelyPairs = avoidLonelyPairs;
        this.exc = null;
    }

    /**
     *
     * @param sequence
     * @param temperature
     */
    public RNAfold(String sequence, double temperature) {
        String command = COMMAND_RNAFOLD; // + RNAFOLD_ARGS;

        InputStreamReader isr = null;
        Process child = null;
        OutputStream out = null;
        BufferedReader br;
        String line;
        String lpCmd = "--noLP";
        String lpTemp = "--temp=";

        this.setTemperature(temperature);

        lpTemp += this.temperature;

        if (!this.avoidLonelyPairs) {
            lpCmd = "";
        }

        try {
            child = new ProcessBuilder(command, "-d2", lpCmd, "--noPS", lpTemp)
                    .start();

            out = child.getOutputStream();
            out.write(sequence.getBytes());
            out.write(13);
            out.write('@');
            out.write(13);
            out.close();

            /*try (*/
            InputStream in = child.getInputStream()/*) {*/;

            isr = new InputStreamReader(in);
            br = new BufferedReader(isr);
            br.readLine();
            br.readLine();
            line = br.readLine();
            this.structure = line.trim()
                    .substring(0, sequence.length());

            Pattern regex = Pattern.compile("(\\d+(?:\\.\\d+)?)");
            Matcher matcher = regex.matcher(line);
            if (matcher.find()) {
                this.mfe = new Double(matcher.group(1));
                this.mfe *= (-1);
            }

            isr.close();
            out.close();
            child.destroy();

        } catch (IOException ex) {
            Logger.getLogger(RNAfold.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("SeqLength: " + sequence.length());
        } finally {
            child = null;
            out = null;
            isr = null;
        }
        /*} catch (Exception ex) {
            err.println("RNAFold error: " + ex.getLocalizedMessage());
            this.exc = ex;
        }*/
    }

    /**
     *
     * @param sequence
     */
    public RNAfold(String sequence) {
        this(sequence, 25.00);
    }

    /**
     *
     * @return the temperature
     */
    public Double getTemperature() {
        return temperature;
    }

    /**
     *
     * @param temperature
     */
    private void setTemperature(Double temperature) {
        this.temperature = temperature;
    }

    public boolean gotError() {
        return exc != null;
    }

    /**
     *
     * @return
     */
    public boolean isAvoidLonelyPairs() {
        return avoidLonelyPairs;
    }

    /**
     *
     * @param avoidLonelyPairs
     */
    public void setAvoidLonelyPairs(boolean avoidLonelyPairs) {
        this.avoidLonelyPairs = avoidLonelyPairs;
    }

    /**
     *
     * @return
     */
    public String getStructure() {
        return structure;
    }

    /**
     *
     * @param structure
     */
    public void setStructure(String structure) {
        this.structure = structure;
    }

    @Override
    /**
     *
     */
    public String toString() {
        return "RnaFoldInterface{" + "structure=" + structure + ", mfe=" + mfe + '}';
    }

    /**
     *
     * @return
     */
    public double getMfe() {
        return mfe;
    }

    /**
     *
     * @param mfe
     */
    public void setMfe(float mfe) {
        this.mfe = mfe;
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {

        String sequence = "UAGAGAUCUCUAUGUAUUUCCC";
        RNAfold test;
        try {
            test = new RNAfold(sequence);
            out.println(test);
        } catch (Exception ex) {
            Logger.getLogger(RNAfold.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
