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

import com.laslo.core.LoopMatcher;
import javax.swing.JButton;
//import javax.swing.JProgressBar;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 *
 * @author David A. Mancilla
 */
class GUISwingWorker extends
        SwingWorker<Integer, Void> {

    private JTextArea textArea;
    private JButton button;
    private LoopMatcher loop;
    //private JProgressBar jp;
    private boolean ok;
    private int progress;

    /**
     * 
     * @param textArea
     * @param button
     * @param loop 
     */
    public GUISwingWorker(JTextArea textArea, JButton button,
            LoopMatcher loop) {
        this.textArea = textArea;
        this.loop = loop;
        this.button = button;
        this.ok = true;
        this.progress = 0;
        //this.jp = jp;
    }

    /**
     * 
     */
    @Override
    protected Integer doInBackground() throws Exception {
        LoopMatcher lm = this.getLoop();
        this.setOk(lm.startReadingFiles());
        //progress = lm.getProgress();
        //setProgress(progress);
        //jp.setValue(progress);
        //System.out.println("task: " + lm.getProgress());
        return 1;
    }

    /**
     * 
     */
    @Override
    protected void done() {
        System.out.flush();
    }

    /**
     * @return the button
     */
    public JButton getButton() {
        return button;
    }

    /**
     * @return the loop
     */
    public LoopMatcher getLoop() {
        return loop;
    }

    /**
     * @return the textArea
     */
    public JTextArea getTextArea() {
        return textArea;
    }

    /**
     * @return the ok
     */
    public boolean isOk() {
        return ok;
    }

    /**
     * @param button the button to set
     */
    public void setButton(JButton button) {
        this.button = button;
    }

    /**
     * @param loop the loop to set
     */
    public void setLoop(LoopMatcher loop) {
        this.loop = loop;
    }

    /**
     * @param ok the ok to set
     */
    public void setOk(boolean ok) {
        this.ok = ok;
    }

    public void setTaskProgress(int progress){
        if(progress < 0) {
            this.progress = 0;
        } else {
            this.progress = progress;
        }
    }
    
    public int getTaskProgress(){
        return this.progress;
    }
    
    /**
     * @param textArea the textArea to set
     */
    public void setTextArea(JTextArea textArea) {
        this.textArea = textArea;
    }
    
}
