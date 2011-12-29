package edu.ucsd.msalign.spec.id;

import edu.ucsd.msalign.spec.align7.preprocess.AlignSpectrum;
import edu.ucsd.msalign.spec.align7.preprocess.AlignSpectrumFactory;
import edu.ucsd.msalign.spec.peak.DeconvPeak;
import edu.ucsd.msalign.spec.sp.Ms;
import java.util.Map;
import java.util.Properties;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Scan;

public class EValueAdapter {

    private static IdCompEValue comp;


    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        init(conf);
        /*
        Map<Integer,Scan> scans = conf.getScans();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Double> evalues = conf.getEvalues();

        int scanId = 3180;
        int proteinId = 4246;
        System.out.println(scanId + " " + proteinId + " " + getBestEValue(scans.get(scanId), proteinId));

        for (int scanId : evalues.keySet()) {
            int proteinId = msAlignResults.get(scanId);
            System.out.println(scanId + " " + proteinId + " " + evalues.get(scanId) + " " + getBestEValue(scans.get(scanId), proteinId));
        }
        */

    }

    public static void init(Configuration conf) throws Exception {
        comp = new IdCompEValue(conf.getProteinDatabaseFile().getCanonicalPath(), conf.getProteinDatabaseFile().getCanonicalPath(), 15);
//        final Properties properties = new Properties();
//        properties.setProperty("databaseFileName", );
//        properties.setProperty("useDecoyDb", "false");
//        properties.setProperty("activation", "CID");
//        properties.setProperty("cysteineProtection", "C57");
//        properties.setProperty("reportNumber", "1");
//        properties.setProperty("shiftNumber", "2");
//        properties.setProperty("errorTolerance", "15");
//        properties.setProperty("eValueThreshold", "0.01");
//        eValueCalculator = new IdEValue();
//        eValueCalculator.init(properties);
    }

    public static double getBestEValue(Scan scan, int proteinId) throws Exception  {
        double best = 9E100;
        Ms<DeconvPeak> deconvSp = scan.getMsDeconvPeaks();
        AlignSpectrum sp[][] = AlignSpectrumFactory.getSpectrum(deconvSp, comp.mng);
        comp.align.setDataSp(scan.getId(), proteinId, sp);
        comp.align.run();
        for (int i = 0; i <= 2; i++) {
            for (int j = 0; j <= 3; j++) {
                double score = comp.getScore(i, j);
                double evalue = comp.getEValue(i, j);
                if (evalue < best && evalue > 0) {
                    best = evalue;
                }
            }
        }

        return best;
    }
}
