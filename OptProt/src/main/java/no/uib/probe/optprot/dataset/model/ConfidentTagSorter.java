/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.dataset.model;

/**
 *
 * @author yfa041
 */
public class ConfidentTagSorter implements Comparable<ConfidentTagSorter> {

    private final Double value;
    private final String title;

    public ConfidentTagSorter(double value, String title) {
        this.value = value;
        this.title = title;
    }

    @Override
    public int compareTo(ConfidentTagSorter o) {
        return this.value.compareTo(o.value);
    }

    public double getValue() {
        return value;
    }

    public String getTitle() {
        return title;
    }

}
