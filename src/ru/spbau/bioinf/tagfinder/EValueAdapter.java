package ru.spbau.bioinf.tagfinder;

import edu.ucsd.msalign.spec.align7.AlignFasta;
import edu.ucsd.msalign.spec.align7.AlignMng;
import edu.ucsd.msalign.spec.align7.AlignResult;
import edu.ucsd.msalign.spec.align7.preprocess.AlignSpectrum;
import edu.ucsd.msalign.spec.align7.preprocess.AlignSpectrumFactory;

public class EValueAdapter {

    private static AlignFasta align;
    private static AlignMng mng;


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
        //comp = new IdCompEValue(conf.getProteinDatabaseFile().getCanonicalPath(), conf.getProteinDatabaseFile().getCanonicalPath(), 15);

        mng = new AlignMng();
        boolean[] newUseSuffix = {true, false, false, false};
        mng.useForSuffixInternal = newUseSuffix;
        mng.internalAllowNShift = new double[1];
        mng.internalAllowNShift[0] = 0.0;
        mng.maxSeqLen = 10000;

        mng.ppo = 15/1000000f;
        mng.nShift = 2;
        mng.resFileName = "residue_extend_carbamidomethylation.xml";
        mng.isRandomSeq = false;
        mng.fragType = AlignMng.CID;

        mng.usePrecMassDelta = false;
        mng.nDelta = 0;
        mng.nTotalDelta = 1;

        align = new AlignFasta(mng);
        align.setDataSeq(conf.getProteinDatabaseFile());

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

    public static double[] getBestPrsm(Scan scan, final int ... proteinIds) throws Exception  {
        double best = 9E100;
        int bestProteinId = -1;
        AlignSpectrum sp[][] = AlignSpectrumFactory.getSpectrum(scan.getMsDeconvPeaks(), mng);
        align.setDataSp(scan.getId(), sp);
        align.updateSeq(proteinIds);
        align.run();
        AlignResult res = align.getBestResult();

        if (res != null) {
            best = res.getEValue();
            bestProteinId = res.getSeqId();
        }

        return new double[] {best, bestProteinId};
    }

}
