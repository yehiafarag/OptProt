/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

import java.util.ArrayList;

/**
 *This class is a java bean class to store optimized search results
 * @author Yehia Mokhtar Farag
 */
public class OptimisedSearchParameters {
    private String digestionParameter;//=DigestionParameters.CleavageParameter.enzyme.name();
    private String  enzymeName="Trypsin";
    private String enzymeSpecificity="specific";
    private int maxMissedCleavage=2;
    private double precursorTolerance;
    private double fragmentTolerance;

    public double getFragmentTolerance() {
        return fragmentTolerance;
    }

    public void setFragmentTolerance(double fragmentTolerance) {
        this.fragmentTolerance = fragmentTolerance;
    }

    public double getPrecursorTolerance() {
        return precursorTolerance;
    }

    public void setPrecursorTolerance(double precursorTolerance) {
        this.precursorTolerance = precursorTolerance;
    }
    private ArrayList<Integer> selectedForwardIons;
    private  ArrayList<Integer> selectedRewindIons ;

    public ArrayList<Integer> getSelectedForwardIons() {
        return selectedForwardIons;
    }

    public void setSelectedForwardIons(ArrayList<Integer> selectedForwardIons) {
        this.selectedForwardIons = selectedForwardIons;
    }

    public ArrayList<Integer> getSelectedRewindIons() {
        return selectedRewindIons;
    }

    public void setSelectedRewindIons(ArrayList<Integer> selectedRewindIons) {
        this.selectedRewindIons = selectedRewindIons;
    }
    

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

    public String getDigestionParameter() {
        return digestionParameter;
    }

    public void setDigestionParameter(String cleavageParameters) {
        this.digestionParameter = cleavageParameters;
    }
}
