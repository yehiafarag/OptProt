/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.configurations;

import java.io.File;

/**
 * Main configurations needed for the sub search utilities
 *
 * @author Yehia Mokhtar Farag
 */
public class Configurations {

    /**
     * The search engine configuration folders.
     */
    public static final String XTANDEM_FOLDER = "D:\\Apps\\searchgui\\resources\\XTandem\\windows\\windows_64bit";
    public static final String MYRIMATCH_FOLDER = "D:\\Apps\\searchgui\\resources\\MyriMatch\\windows\\windows_64bit";
    public static final String NOVOR_FOLDER = "D:\\Apps\\searchgui\\resources\\Novor";
    public static final String COMET_FOLDER = "D:\\Apps\\searchgui\\resources\\Comet\\windows";
    public static final String DIRECTAG_FOLDER = "D:\\Apps\\searchgui\\resources\\DirecTag\\windows\\windows_64bit";
    /**
     * The configurations folder.
     */
//    public static final String CONFIG_FOLDER = null;
    /**
     * The resources folder.
     */
//    public static final String RESOURCE_FOLDER = "D:\\Apps\\OptProt\\resources\\";
    public static final String DEFAULT_RESULT_NAME = "optsearch_results";

    /**
     * The output folder.
     */
    public static final String OUTPUT_FOLDER_PATH = "D:\\Apps\\OptProt\\data\\output";
    /**
     * The default search param file.
     */
    public static final String DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE = "D:\\Apps\\OptProt\\data\\default_optprot_search_settings.par";

    /**
     * The active search param file.
     */
    public static File ACTIVE_SEARCH_SETTINGS_FILE;

    public static String Dataset_Id;

    public static final String EXTRACT_MS_TYPE = "TA";//TA  WF
//    public static int EXTRACT_MAX_MS_SIZE = 1500;
//    public static int EXTRACT_MIN_MS_SIZE = 1000;
//    public static int MIN_TAG_SIZE = 1000;
    public static int REFINED_MS_SIZE =2000;

    public static String get_current_file_fingerprent() {
        return "_" + EXTRACT_MS_TYPE;
    }

}
