package com.tools.joewandy;

public class Protein implements Sequence {

	private String sequenceString;
	private String id;

	public Protein(String id, String sequenceString) {
		this(sequenceString);
		this.id = id;
	}

	public Protein(String sequenceString) {
		this.sequenceString = sequenceString.replaceAll("\\s", ""); //$NON-NLS-1$ //$NON-NLS-2$
	}

	@Override
	public String getSequenceString() {
		return this.sequenceString;
	}

	@Override
	public void setSequenceString(String sequenceString) {
		this.sequenceString = sequenceString;
	}

	@Override
	public String getId() {
		return this.id;
	}

	@Override
	public void setId(String id) {
		this.id = id;
	}

	@Override
	public String toString() {
		return "Dna [id=" + this.id + ", sequenceString=" + this.sequenceString + "]"; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
	}

}
