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
 * @author David A. Mancilla
 *
 */
public class BioMartFastaID extends SourceFile {

    /**
     *
     */
    protected static int nCols = 6;

    /**
     *
     */
    protected static char fs = '|';

    /**
     *
     */
    protected static String columns[] = new String[nCols];

    /**
     *
     * @return
     */
    public static int getNumColumns() {
        return columns.length;
    }

    /**
     *
     * @param cols
     */
    public static void setColumns(String cols[]) {
        columns = cols;
    }

    /**
     *
     * @return
     */
    public static String getHeader() {
        String header = "";

        for (int i = 0; i < nCols; i++) {
            header += "Column" + (i + 1) + ROW_DELIMITER;
        }

        return header;
    }

    /**
     *
     * @param idSequence
     */
    public void setBioMartTags(String idSequence) {

        String auxSeq = "";

        if (idSequence.contains("|")) {
            auxSeq = idSequence.replace(fs, '@');
        } else {
            auxSeq = idSequence.replace('=', '@');
        }

        auxSeq = auxSeq.replace(';', '-');
        String[] cols = auxSeq.split("@");

        for (int i = 0; i < cols.length && i < 6; i++) {
            columns[i] = cols[i];
        }
    }

    /**
     *
     * @return
     */
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
