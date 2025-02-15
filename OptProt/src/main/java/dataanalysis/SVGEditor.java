/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package dataanalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yfa041
 */
public class SVGEditor {

    private File resultsFolder;

    public void handelVennFiles() {
        File srcFolder = new File("C:\\Users\\yfa041\\OneDrive - University of Bergen\\manuscript 2024\\Analysis");
        resultsFolder = new File(srcFolder, "svgResults");
        File ratioFile = new File(srcFolder, "Ratios.csv");

        System.out.println("Source folder exist " + srcFolder.exists() + "  and ratio file " + ratioFile.exists());
        Map<String, String[]> datasetIntersectionMap = new LinkedHashMap<>();
        try {
            FileReader reader = new FileReader(ratioFile);
            BufferedReader bufferedReader = new BufferedReader(reader);

            String line;

            while ((line = bufferedReader.readLine()) != null) {
                if (line.startsWith("pxd")) {
                    String[] dsInfor = line.split(";");
                    // sage(full) intersection	gained	lost	PSMs Overall gained - (pride) intersection  gained	lost	PSMs Overall gained   *** extan (full) intersection	gained	lost	PSMs Overall gained - (pride) intersection  gained	lost	PSMs Overall gained timeSageSub timeSagefull timeXtanSub timeXtanFull

                    //                  0           1     2                 3                       4           5          6            7                                   8               9       10          11                      12              13          14          15              16          17          18          19      
                    datasetIntersectionMap.put(dsInfor[0], new String[]{(dsInfor[3]), (dsInfor[4]), (dsInfor[5]), (dsInfor[6]), (dsInfor[8]), (dsInfor[9]), (dsInfor[10]), (dsInfor[11]), (dsInfor[14]), (dsInfor[15]), (dsInfor[16]), (dsInfor[17]), (dsInfor[19]), dsInfor[20], dsInfor[21], dsInfor[22], dsInfor[23], dsInfor[24], dsInfor[25], dsInfor[26]});
                }
            }
            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

        for (String datasetName : datasetIntersectionMap.keySet()) {
            File datasetFolder = new File(srcFolder, datasetName);
            for (File graphFile : datasetFolder.listFiles()) {
                if (graphFile.getName().endsWith(".svg")) {
//                    System.out.println("Svg file exist in " + datasetName + "  " + graphFile.getName());

// sage(full) intersection	gained	lost	PSMs Overall gained - (pride) intersection  gained	lost	PSMs Overall gained   *** extan (full) intersection	gained	lost	PSMs Overall gained - (pride) intersection  gained	lost	PSMs Overall gained timeSageSub timeSagefull timeXtanSub timeXtanFull
//                  0           1     2                 3                       4           5          6            7                                   8               9       10          11                      12              13          14          15              16          17          18          19      
                    if (graphFile.getName().contains("full") && graphFile.getName().contains("sage")) {
                        editExportGraphFile(datasetName, graphFile, datasetIntersectionMap.get(datasetName)[0], datasetIntersectionMap.get(datasetName)[1], datasetIntersectionMap.get(datasetName)[2], datasetIntersectionMap.get(datasetName)[3], datasetIntersectionMap.get(datasetName)[16], datasetIntersectionMap.get(datasetName)[17]);
                    } else if (graphFile.getName().contains("pride") && graphFile.getName().contains("sage")) {
                        editExportGraphFile(datasetName, graphFile, datasetIntersectionMap.get(datasetName)[4], datasetIntersectionMap.get(datasetName)[5], datasetIntersectionMap.get(datasetName)[6], datasetIntersectionMap.get(datasetName)[7], "", "");
                    } else if (graphFile.getName().contains("full") && graphFile.getName().contains("xtand")) {
                        editExportGraphFile(datasetName, graphFile, datasetIntersectionMap.get(datasetName)[8], datasetIntersectionMap.get(datasetName)[9], datasetIntersectionMap.get(datasetName)[10], datasetIntersectionMap.get(datasetName)[11], datasetIntersectionMap.get(datasetName)[18], datasetIntersectionMap.get(datasetName)[19]);
                    } else if (graphFile.getName().contains("pride") && graphFile.getName().contains("xtand")) {
                        editExportGraphFile(datasetName, graphFile, datasetIntersectionMap.get(datasetName)[12], datasetIntersectionMap.get(datasetName)[13], datasetIntersectionMap.get(datasetName)[14], datasetIntersectionMap.get(datasetName)[15], "", "");
                    } else {
                        System.out.println("--------error new file name error------" + graphFile.getName());
                    }
                }
                try {
                    Thread.sleep(100);
                } catch (InterruptedException ex) {
                    Logger.getLogger(SVGEditor.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

        }
    }

    private void editExportGraphFile(String dataset_id, File svgInputFile, String intersection, String subGained, String lost, String overall, String subTime, String fullTime) {
        StringBuilder sb = new StringBuilder();
        String fullLine = "";
        try {
            try (FileReader reader = new FileReader(svgInputFile)) {
                BufferedReader bufferedReader = new BufferedReader(reader);                
                String line;                
                while ((line = bufferedReader.readLine()) != null) {
                    fullLine += line;
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        String path = ""; 
        String fullPride = "Full Data Mode";
        if (fullLine.contains("pri")) {
            fullPride = "PRIDE Mode";
        }
        fullLine = fullLine.replace("<svg width=\"400\" height=\"400\"", "<svg width=\"500\" height=\"450\"").replace("fill-opacity: 0.3; fill: rgb(31, 119, 180)", "fill-opacity: 0.5; fill: rgb(0, 176, 240)").replace("rgb(255, 127, 14)", "rgb(255, 121, 0)").replace("pri", intersection + "");
        String[] circls = fullLine.split("<g><circle");
        sb.append(circls[0]);
        sb.append("<g><circle ");
        //circle 1 subset 
        if (circls[1].contains("fill-opacity: 0.3;")) {
            String first = circls[2];
            String second = circls[1];
            path = "<path " + first.split("<path")[1];
            first = first.split("<path")[0];
            circls[1] = first;
            circls[2] = second;
        }
///     edit location
        String editCircleI = circls[1]; 
        String cx1 = editCircleI.split("cx=")[1].split("cy=")[0];
        String updatedCx1 = "\"" + (50.0 + Double.parseDouble(cx1.replace("\"", "")) + "\" ");
        editCircleI = editCircleI.replace(cx1, updatedCx1);
        //edit label 
        String subLabelText = "<text" + editCircleI.split("<text")[1].split("</text>")[0] + "</text>"; 
        String updatedSubLabelText = "<text dy=\".35em\" x=\"55\" y=\"25\" text-anchor=\"middle\"  style=\"fill: rgb(0, 176, 240); font-family: Trebuchet MS, Tahoma, sans-serif; font-size: 10pt;\"> <tspan x=\"55\" y=\"25\" dy=\"0.35em\">Subset Mode</tspan><tspan x=\"55\" y=\"45\" dy=\"0.35em\">Gained PSMs #</tspan>   <tspan x=\"55\" y=\"65\" dy=\"0.35em\">(+) " + subGained + "</tspan> <tspan x=\"50\" y=\"85\" dy=\"0.35em\">" + subTime + "</tspan></text>";
        editCircleI = editCircleI.replace(subLabelText, updatedSubLabelText);
        sb.append(editCircleI);
        sb.append("<g><circle");
        String editCircleII = circls[2];
        String cx2 = editCircleII.split("cx=")[1].split("cy=")[0];
        String updatedCx2 = "\"" + (50.0 + Double.parseDouble(cx2.replace("\"", "")) + "\" ");
        editCircleII = editCircleII.replace(cx2, updatedCx2);
        String otherLabelText = "<text" + editCircleII.split("<text")[1].split("</text>")[0] + "</text>";
        String updatedOtherLabelText = " <text dy=\".35em\" x=\"445\" y=\"20\" text-anchor=\"middle\" style=\"fill: rgb(255, 121, 0); font-family: Trebuchet MS, Tahoma, sans-serif; font-size: 10pt;\"> <tspan x=\"445\" y=\"25\" dy=\"0.35em\">" + fullPride + "</tspan><tspan x=\"445\" y=\"45\" dy=\"0.35em\">Lost PSMs #</tspan> <tspan x=\"445\" y=\"65\" dy=\"0.35em\">(-) " + lost + "</tspan> <tspan x=\"445\" y=\"85\" dy=\"0.35em\">" + fullTime + "</tspan></text>";
        editCircleII = editCircleII.replace(otherLabelText, updatedOtherLabelText);//                  
        sb.append(editCircleII);
        sb.append(path);
        String extraText = "<text font-weight=\"bold\" dy=\".35em\" x=\"250\" y=\"180\" text-anchor=\"middle\" style=\"fill: rgb(255, 255, 255); font-family: Trebuchet MS, Tahoma, sans-serif; font-size: 12pt;\"><tspan x=\"250\" y=\"190\" dy=\"0.35em\">Intersection</tspan><tspan x=\"250\" y=\"215\" dy=\"0.35em\">PSMs # " + intersection + "</tspan></text>" + "<text  width=\"500\" dy=\".35em\" x=\"250\" y=\"400\" text-anchor=\"middle\" style=\"fill: rgb(80,80,80); font-family: Trebuchet MS, Tahoma, sans-serif; font-size: 12pt;\"><tspan x=\"250\" y=\"410\" dy=\"0.35em\">" + dataset_id.toUpperCase() + "</tspan><tspan x=\"250\" y=\"430\" dy=\"0.35em\">Overall gained PSMs # " + overall + "</tspan></text>" + "</svg>";
        String svgAsText = sb.toString().replace("</svg>", extraText);
        try {
            File outputFile = new File(resultsFolder, dataset_id + "_" + svgInputFile.getName());
            if (outputFile.exists()) {
                outputFile.delete();
            }
            outputFile.createNewFile();
            FileWriter writer = new FileWriter(outputFile, false);
            try (BufferedWriter bufferedWriter = new BufferedWriter(writer)) {
                bufferedWriter.write(svgAsText);
                bufferedWriter.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        SVGEditor svgEditor = new SVGEditor();
        svgEditor.handelVennFiles();
    }
}
