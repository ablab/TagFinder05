package ru.spbau.bioinf.tagfinder;


import edu.ucsd.msalign.align.MsAlign;
import edu.ucsd.msalign.align.PropertyUtil;
import edu.ucsd.msalign.align.prsm.PrSM;
import java.util.List;
import java.util.Properties;


public class EValueAdapter {

    private static MsAlign align;


    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        init(conf);
    }

    public static void init(Configuration conf) throws Exception {
        //comp = new IdCompEValue(conf.getProteinDatabaseFile().getCanonicalPath(), conf.getProteinDatabaseFile().getCanonicalPath(), 15);
        conf.getScans();
        Properties prop = PropertyUtil.getDefaultProperties();
        prop.put("databaseFileName", conf.getProteinDatabaseFile().getAbsolutePath());
        prop.put("spectrumFileName", conf.getMsalignFile().getAbsolutePath());
        prop.put("activation", "CID");
        prop.put("cysteineProtection", "C57");
        align = new MsAlign(prop);
        align.readSpectra();
    }

    public static double[] getBestPrsm(Scan scan, final int ... proteinIds) throws Exception  {
        double best = 9E100;
        int bestProteinId = -1;
        align.updateSeqs(proteinIds);
        List<PrSM[]> prsms = align.getPrsms(scan.getId());

        for (PrSM[] prsm : prsms) {
            for (PrSM p : prsm) {
                if (p.getEValue() < best) {
                    best = p.getEValue();
                    bestProteinId = p.getSeq().getId();
                }                
            }
        } 

        return new double[] {best, bestProteinId};
    }

}
