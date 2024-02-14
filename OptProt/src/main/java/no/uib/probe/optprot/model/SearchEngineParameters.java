/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

import java.io.File;

/**
 *
 * @author yfa041
 */
public class SearchEngineParameters {
     private boolean runOmssa = false;
            private boolean runXTandem = true;
            private boolean runMsgf = false;
            private boolean runMsAmanda = false;
            private boolean runMyriMatch = false;
            private boolean runComet = false;
            private boolean runTide = false;
            private boolean runAndromeda = false;
            private boolean runMetaMorpheus = false;
            private boolean runSage = false;
            private boolean runNovor = false;
            private boolean runDirecTag = false;
            private File omssaFolder = null;
            private File xTandemFolder =null;
            private File msgfFolder = null;
            private File msAmandaFolder = null;
            private File myriMatchFolder = null;
            private File cometFolder = null;
            private File tideFolder = null;
            private File tideIndexLocation = null;
            private File andromedaFolder = null;
            private File metaMorpheusFolder = null;
            private File sageFolder = new File("D:\\Apps\\OptProt\\data\\saga");
            private File novorFolder = null;
            private File direcTagFolder = null;
            private File makeblastdbFolder = null;

    public boolean isRunOmssa() {
        return runOmssa;
    }

    public void setRunOmssa(boolean runOmssa) {
        this.runOmssa = runOmssa;
    }

    public boolean isRunXTandem() {
        return runXTandem;
    }

    public void setRunXTandem(boolean runXTandem) {
        this.runXTandem = runXTandem;
    }

    public boolean isRunMsgf() {
        return runMsgf;
    }

    public void setRunMsgf(boolean runMsgf) {
        this.runMsgf = runMsgf;
    }

    public boolean isRunMsAmanda() {
        return runMsAmanda;
    }

    public void setRunMsAmanda(boolean runMsAmanda) {
        this.runMsAmanda = runMsAmanda;
    }

    public boolean isRunMyriMatch() {
        return runMyriMatch;
    }

    public void setRunMyriMatch(boolean runMyriMatch) {
        this.runMyriMatch = runMyriMatch;
    }

    public boolean isRunComet() {
        return runComet;
    }

    public void setRunComet(boolean runComet) {
        this.runComet = runComet;
    }

    public boolean isRunTide() {
        return runTide;
    }

    public void setRunTide(boolean runTide) {
        this.runTide = runTide;
    }

    public boolean isRunAndromeda() {
        return runAndromeda;
    }

    public void setRunAndromeda(boolean runAndromeda) {
        this.runAndromeda = runAndromeda;
    }

    public boolean isRunMetaMorpheus() {
        return runMetaMorpheus;
    }

    public void setRunMetaMorpheus(boolean runMetaMorpheus) {
        this.runMetaMorpheus = runMetaMorpheus;
    }

    public boolean isRunSage() {
        return runSage;
    }

    public void setRunSage(boolean runSage) {
        this.runSage = runSage;
    }

    public boolean isRunNovor() {
        return runNovor;
    }

    public void setRunNovor(boolean runNovor) {
        this.runNovor = runNovor;
    }

    public boolean isRunDirecTag() {
        return runDirecTag;
    }

    public void setRunDirecTag(boolean runDirecTag) {
        this.runDirecTag = runDirecTag;
    }

    public File getOmssaFolder() {
        return omssaFolder;
    }

    public void setOmssaFolder(File omssaFolder) {
        this.omssaFolder = omssaFolder;
    }

    public File getxTandemFolder() {
        return xTandemFolder;
    }

    public void setxTandemFolder(File xTandemFolder) {
        this.xTandemFolder = xTandemFolder;
    }

    public File getMsgfFolder() {
        return msgfFolder;
    }

    public void setMsgfFolder(File msgfFolder) {
        this.msgfFolder = msgfFolder;
    }

    public File getMsAmandaFolder() {
        return msAmandaFolder;
    }

    public void setMsAmandaFolder(File msAmandaFolder) {
        this.msAmandaFolder = msAmandaFolder;
    }

    public File getMyriMatchFolder() {
        return myriMatchFolder;
    }

    public void setMyriMatchFolder(File myriMatchFolder) {
        this.myriMatchFolder = myriMatchFolder;
    }

    public File getCometFolder() {
        return cometFolder;
    }

    public void setCometFolder(File cometFolder) {
        this.cometFolder = cometFolder;
    }

    public File getTideFolder() {
        return tideFolder;
    }

    public void setTideFolder(File tideFolder) {
        this.tideFolder = tideFolder;
    }

    public File getTideIndexLocation() {
        return tideIndexLocation;
    }

    public void setTideIndexLocation(File tideIndexLocation) {
        this.tideIndexLocation = tideIndexLocation;
    }

    public File getAndromedaFolder() {
        return andromedaFolder;
    }

    public void setAndromedaFolder(File andromedaFolder) {
        this.andromedaFolder = andromedaFolder;
    }

    public File getMetaMorpheusFolder() {
        return metaMorpheusFolder;
    }

    public void setMetaMorpheusFolder(File metaMorpheusFolder) {
        this.metaMorpheusFolder = metaMorpheusFolder;
    }

    public File getSageFolder() {
        return sageFolder;
    }

    public void setSageFolder(File sageFolder) {
        this.sageFolder = sageFolder;
    }

    public File getNovorFolder() {
        return novorFolder;
    }

    public void setNovorFolder(File novorFolder) {
        this.novorFolder = novorFolder;
    }

    public File getDirecTagFolder() {
        return direcTagFolder;
    }

    public void setDirecTagFolder(File direcTagFolder) {
        this.direcTagFolder = direcTagFolder;
    }

    public File getMakeblastdbFolder() {
        return makeblastdbFolder;
    }

    public void setMakeblastdbFolder(File makeblastdbFolder) {
        this.makeblastdbFolder = makeblastdbFolder;
    }
}