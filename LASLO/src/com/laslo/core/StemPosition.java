/**
 * 
 */
package com.laslo.core;

/**
 * @author David
 *
 */
public class StemPosition {
	String id;
	int position;
	int range;
	StemPosition(String id, int pos){
		this.id = id;
		this.position = pos;
		this.range = 12;
	}

	public int compareTo(StemPosition arg0) {
		// TODO Auto-generated method stub
		int ret = -1, val;
		
		if(arg0.id.equals(this.id)){
			val = Math.abs(this.position - arg0.position);
			
			if(val <= range)
				ret = 0;
		}
		
		return ret;
	}

	public int compare(StemPosition arg0, StemPosition arg1) {
		// TODO Auto-generated method stub
		return arg0.compareTo(arg1);
	}

}
