/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import java.io.File;
import no.uib.probe.optprot.model.OptProtSearchParameters;

/**
 *
 * @author yfa041
 */
public class SearchOptimizerUtilities {
    public static  OptProtSearchParameters initSearchOptimizerParameters() {
        OptProtSearchParameters searchOptimizerParameters = new OptProtSearchParameters();
        searchOptimizerParameters.setxTandemFolder(new File("D:\\Apps\\searchgui\\resources\\XTandem\\windows\\windows_64bit"));
        searchOptimizerParameters.setNovorFolder(new File("D:\\Apps\\searchgui\\resources\\Novor"));
        searchOptimizerParameters.setCometFolder(new File("D:\\Apps\\searchgui\\resources\\Comet\\windows"));
        searchOptimizerParameters.setDirecTagFolder(new File("D:\\Apps\\searchgui\\resources\\DirecTag\\windows\\windows_64bit"));
        return searchOptimizerParameters;
    }
}
