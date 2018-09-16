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
public class GenericID extends SourceFile {

    final static String HEADER = "SequenceID" + ROW_DELIMITER;

    /**
     *
     */
    public GenericID() {
        this.geneID = "";
    }

    @Override
    public String toString() {
        return "GenericID [sequenceID=" + geneID + "]";
    }

    /**
     *
     * @return
     */
    public static String getHeader() {
        return GenericID.HEADER;
    }

    /**
     *
     * @param id
     */
    public void setGenericTags(String id) {
        // TODO Auto-generated method stub
        geneID = id.replace(';', '-');
    }

    @Override
    public String toRowCSV() {
        return geneID + ROW_DELIMITER;
    }
}
