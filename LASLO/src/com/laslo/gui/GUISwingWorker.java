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
package com.laslo.gui;

import com.laslo.core.LoopCatcher;
import javax.swing.JButton;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 *
 * @author David A. Mancilla
 */
class GUISwingWorker extends
        SwingWorker<Integer, Void> {

    private final JTextArea textArea;
    private final JButton button;
    private final LoopCatcher loop;
    private boolean ok;

    /**
     * 
     * @param textArea
     * @param button
     * @param loop 
     */
    public GUISwingWorker(JTextArea textArea, JButton button,
            LoopCatcher loop) {
        this.textArea = textArea;
        this.loop = loop;
        this.button = button;
        this.ok = true;
    }

    /**
     * 
     */
    @Override
    protected Integer doInBackground() throws Exception {
        this.ok = this.loop.beginSearch();
        return 1;
    }

    /**
     * 
     */
    @Override
    protected void done() {
        System.out.flush();
        //button.setEnabled(false);
    }

    /*@Override
     protected void process(List<Object[]> chunks) {
         for (Object[] row : chunks) {
             tableModel.addRow(row);
         }
     }*/
}
