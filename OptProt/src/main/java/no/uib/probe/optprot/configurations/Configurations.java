/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.configurations;

/**
 * Main configurations needed for the sub search utilities
 *
 * @author Yehia Mokhtar Farag
 */
public class Configurations {

    /**
     * The configurations folder.
     */
    public static final String XTANDEM_FOLDER = "D:\\Apps\\searchgui\\resources\\XTandem\\windows\\windows_64bit";
    /**
     * The configurations folder.
     */
//    public static final String CONFIG_FOLDER = null;
    /**
     * The resources folder.
     */
    public static final String RESOURCE_FOLDER = "D:\\Apps\\OptProt\\resources\\";
    public static final String DEFAULT_RESULT_NAME = "optsearch_results";

    /**
     * The output folder.
     */
    public static final String OUTPUT_FOLDER_PATH = "D:\\Apps\\OptProt\\data\\output";

    public static String Dataset_Id;

    public static final String EXTRACT_MS_TYPE = "TA";//TA  WF
    public static int EXTRACT_MAX_MS_SIZE = 3500;
    public static int EXTRACT_MIN_MS_SIZE = 3000;
    public static int MIN_TAG_SIZE = 1000;

    public static String get_current_file_fingerprent() {
        return "_" + EXTRACT_MS_TYPE;
    }

}
