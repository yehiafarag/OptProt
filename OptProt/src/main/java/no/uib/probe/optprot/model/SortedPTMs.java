/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

/**
 *
 * @author yfa041
 */
public class SortedPTMs implements Comparable<SortedPTMs> {

    private final String name;
    private final double score;

    public SortedPTMs(String name, double score) {
        this.name = name;
        this.score = score;
    }

    @Override
    public int compareTo(SortedPTMs o) {
        return Double.compare(score, o.score);
    }

    @Override
    public String toString() {
        return name+"  "+score;
    }

    public String getName() {
        return name;
    }

    public double getScore() {
        return score;
    }

}
