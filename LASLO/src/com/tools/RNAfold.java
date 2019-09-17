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
public class RNAFold {

    private final static String COMMAND_RNAFold = "./ext/RNAFold.exe";
    //private final static String RNAFold_ARGS = " -d2 --noLP --noPS";
    private String structure;
    private double mfe;
    private Double temperature;
    private boolean avoidLonelyPairs;
    private Exception exc;

    /**
     *
     */
    public RNAFold() {
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
    public RNAFold(double temp, boolean avoidLonelyPairs) {
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
    public RNAFold(String sequence, double temperature, boolean avoidLonelyPairs) {
        String command = COMMAND_RNAFold; // + RNAFold_ARGS;

        InputStreamReader isr = null;
        Process child = null;
        OutputStream outstr = null;
        BufferedReader br;
        String line;
        String lpCmd = "--noLP";
        String lpTemp = "--temp=";

        this.setTemperature(temperature);

        lpTemp += this.temperature;

        this.avoidLonelyPairs = avoidLonelyPairs;

        if (!this.avoidLonelyPairs) {
            lpCmd = "";
        }

        try {
            child = new ProcessBuilder(command, "-d2", lpCmd, "--noPS", lpTemp)
                    .start();

            outstr = child.getOutputStream();
            outstr.write(sequence.getBytes());
            outstr.write(13);
            outstr.write('@');
            outstr.write(13);
            outstr.close();

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
            outstr.close();
            child.destroy();

        } catch (IOException ex) {
            out.println("--Error executing RNAFold--");
            //Logger.getLogger(RNAFold.class.getName()).log(Level.SEVERE, null, ex);
            out.println("SeqLength: " + sequence.length());
        } finally {
            child = null;
            outstr = null;
            isr = null;
        }
    }

    /**
     *
     * @param sequence
     */
    public RNAFold(String sequence) {
        this(sequence, 25.00, true);
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
        return "RNAFoldInterface{" + "structure=" + structure + ", mfe=" + mfe + '}';
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
        RNAFold test;
        try {
            test = new RNAFold(sequence);
            out.println(test);
        } catch (Exception ex) {
            Logger.getLogger(RNAFold.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}