package no.uib.probe.optprot.util;

import com.compomics.util.experiment.biology.proteins.Protein;
import com.compomics.util.experiment.io.biology.protein.Header;
import com.compomics.util.experiment.io.biology.protein.iterators.FastaIterator;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileWriter;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.io.IoUtil;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yfa041
 */
public class SpectraFileUtilities {

    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param mgfFile The spectrum file to use to get the spectra.
     * @param destinationFile The file where to write.
     * @param maxSpectNum number of spectra in the final file
     */
    public static void writeSubSetEveryNMgfFile(
            File mgfFile,
            File destinationFile,
            double maxSpectNum
    ) {
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(mgfFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(mgfFile, new OptProtWaitingHandler());
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            int n = (int) Math.round(spectrumTitles.length / maxSpectNum);
            if (spectrumTitles == null) {
                throw new IllegalArgumentException(
                        fileNameWithoutExtension + " not loaded."
                );
            }
            MgfFileWriter writer = new MgfFileWriter(destinationFile);
            for (int i = 0; i < spectrumTitles.length;) {
                String spectrumTitle = spectrumTitles[i];
                Spectrum spectrum = msFileHandler.getSpectrum(fileNameWithoutExtension, spectrumTitle);
                writer.writeSpectrum(spectrumTitle, spectrum);
                i += n;
            }

            writer.close();
        } catch (IOException ex) {
            Logger.getLogger(SpectraFileUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param mgfFile The spectrum file to use to get the spectra.
     * @param destinationFile The file where to write.
     *
     * @param maxSpectNum number of spectra in the final file
     */
    public static void writeSubSetTargtedFileAreaMgfFile(
            File mgfFile,
            File destinationFile,
            double maxSpectNum
    ) {
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(mgfFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(mgfFile, new OptProtWaitingHandler());
            
            
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            int targetedAreaStep = spectrumTitles.length / 4;
            int n = (int) Math.round(targetedAreaStep * 2 / maxSpectNum);
            
            if (spectrumTitles == null) {
                throw new IllegalArgumentException(
                        fileNameWithoutExtension + " not loaded."
                );
            }
            MgfFileWriter writer = new MgfFileWriter(destinationFile);
            for (int i = targetedAreaStep; i < (spectrumTitles.length - targetedAreaStep);) {
                String spectrumTitle = spectrumTitles[i];
                Spectrum spectrum = msFileHandler.getSpectrum(fileNameWithoutExtension, spectrumTitle);
                writer.writeSpectrum(spectrumTitle, spectrum);
                i += n;
            }

            writer.close();
            System.out.println("at n: " + n + "  " + targetedAreaStep * 2 + "  " + spectrumTitles.length);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Appends decoy sequences to the provided FASTA file.
     *
     * @param fastaIn the FASTA file to read
     * @param fastaOut the FASTA file to write displaying progress
     * @param micro
     */
    public static void createSubFastaFile(
            File fastaIn,
            File fastaOut,
            boolean micro
    ) {

        try {
            int total = 0;
            int left = 0;
            int specIndex = 0;
            FastaIterator fastaIterator = new FastaIterator(fastaIn);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {
                Protein protein;
                while ((protein = fastaIterator.getNextProtein()) != null) {
                    Header header = fastaIterator.getLastHeader();
                    total++;
                    if (header.getProteinEvidence() != 1 || protein.getLength() <= 400) { //&& !header.asGenericHeader().contains("SV=2")
                        continue;
                    }
                    specIndex++;
                    if (micro) {
                        if (specIndex == 3) {
                            specIndex = 0;
                            continue;
                        }
                    }
                    left++;
                    String accession = protein.getAccession();
                    String sequence = protein.getSequence();

                    String rawHeader = header.getRawHeader();

                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);

                    bw.newLine();
//                    bw.newLine();
//
//                    int accessionEndIndex = rawHeader.indexOf(accession) + accession.length();
//
//                    String part0 = rawHeader.substring(0, rawHeader.indexOf(accession));
//                    String part1 = rawHeader.substring(rawHeader.indexOf(accession), accessionEndIndex);
//                    String part2 = rawHeader.substring(accessionEndIndex);
//
//                    bw.write(part0);
//
//                    bw.write(part1);
//
//                    bw.write(part2);
//
//                    bw.newLine();
//
//                    char[] sequenceAsArray = protein.getSequence().toCharArray();
//
//                    for (int i = sequenceAsArray.length - 1; i >= 0; i--) {
//
//                        char aa = sequenceAsArray[i];
//
//                        bw.write(aa);
//
//                    }

//                    bw.newLine();
//                    bw.newLine();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            System.out.println("total: " + total + "  left " + left);
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }

    public static void createSubFastaFile(
            File fastaIn,
            File fastaOut,
            Set<String> sequences
    ) {

        try {
            int total = 0;
            int left = 0;
            int specIndex = 0;
            FastaIterator fastaIterator = new FastaIterator(fastaIn);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {
                Protein protein;
                while ((protein = fastaIterator.getNextProtein()) != null) {
                    Header header = fastaIterator.getLastHeader();
                    total++;
                    if (header.getProteinEvidence() != 1) { //&& !header.asGenericHeader().contains("SV=2")|| protein.getLength() <= 400
                        continue;
                    }
                    specIndex++;
                    boolean escape = true;
                    for (String seq : sequences) {
                        if (protein.getSequence().contains(seq)) {
                            escape = false;
                            break;
                        }
//                        System.out.println("at prot sec "+seq+"  --- "+protein.getSequence()+" ");

                    }
                    if (escape) {
                        continue;
                    }

                    left++;
                    String accession = protein.getAccession();
                    String sequence = protein.getSequence();

                    String rawHeader = header.getRawHeader();

                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);

                    bw.newLine();
//                    bw.newLine();
//
//                    int accessionEndIndex = rawHeader.indexOf(accession) + accession.length();
//
//                    String part0 = rawHeader.substring(0, rawHeader.indexOf(accession));
//                    String part1 = rawHeader.substring(rawHeader.indexOf(accession), accessionEndIndex);
//                    String part2 = rawHeader.substring(accessionEndIndex);
//
//                    bw.write(part0);
//
//                    bw.write(part1);
//
//                    bw.write(part2);
//
//                    bw.newLine();
//
//                    char[] sequenceAsArray = protein.getSequence().toCharArray();
//
//                    for (int i = sequenceAsArray.length - 1; i >= 0; i--) {
//
//                        char aa = sequenceAsArray[i];
//
//                        bw.write(aa);
//
//                    }

//                    bw.newLine();
//                    bw.newLine();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            System.out.println("total: " + total + "  left " + left);
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }

    public static Set<String> getSequences(File novorOutputFile) {
        String line;
        String splitBy = ",";
        Set<String> sequences = new HashSet<>();
        try {
//parsing a CSV file into BufferedReader class constructor  
            BufferedReader br = new BufferedReader(new FileReader(novorOutputFile));
            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                if (line.startsWith("# id,")) {
                    break;
                }
            }

            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                String[] values = line.split(splitBy);    // use comma as separator  
                sequences.add(values[9].replace("(0)", "").trim());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequences;
    }

}
