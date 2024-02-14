/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

/**
 *This class is a java bean class to store optimized search results
 * @author Yehia Mokhtar Farag
 */
public class OptimisedSearchParameters {
    private String cleavageParameters;// =DigestionParameters.CleavageParameter.enzyme.name();
    private String  enzymeName;//="Trypsin";
    private String enzymeSpecificity;//="specific";
    private int maxMissedCleavage=-1;

    public int getMaxMissedCleavage() {
        return maxMissedCleavage;
    }

    public void setMaxMissedCleavage(int maxMissedCleavage) {
        this.maxMissedCleavage = maxMissedCleavage;
    }

    public String getEnzymeSpecificity() {
        return enzymeSpecificity;
    }

    public void setEnzymeSpecificity(String enzymeSpecificity) {
        this.enzymeSpecificity = enzymeSpecificity;
    }

    public String getEnzymeName() {
        return enzymeName;
    }

    public void setEnzymeName(String enzymeName) {
        this.enzymeName = enzymeName;
    }

    public String getCleavageParameters() {
        return cleavageParameters;
    }

    public void setCleavageParameters(String cleavageParameters) {
        this.cleavageParameters = cleavageParameters;
    }
}
