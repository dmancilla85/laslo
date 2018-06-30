/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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

     public GUISwingWorker(JTextArea textArea, JButton button,
             LoopCatcher loop) {
         this.textArea = textArea;
         this.loop = loop;
         this.button = button;
         this.ok = true;
     }

     @Override
     protected Integer doInBackground() throws Exception {
        this.ok = this.loop.beginResearch();
        return 1;
     }

     @Override
     protected void done(){
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