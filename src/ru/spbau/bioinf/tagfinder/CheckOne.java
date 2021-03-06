package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CheckOne {
    private static NumberFormat df = NumberFormat.getInstance();
    private static List<Protein> proteins;
    private static Configuration conf;

    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        //Logger.getLogger(AlignFasta.class).setLevel(Level.TRACE);
        //Logger.getLogger(AlignEValue.class).setLevel(Level.TRACE);
        //Logger.getLogger(CompPValue.class).setLevel(Level.TRACE);
        //Logger.getLogger(AlignSpectrum.class).setLevel(Level.TRACE);
        long start = System.currentTimeMillis();                
        conf = new Configuration(args);
        EValueAdapter.init(conf, "TARGET");
        proteins = conf.getProteins();
        Tag4Finder.proteins = proteins;
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        //keys.add(3695); keys.add(2239); keys.add(1718); keys.add(1264); keys.add(1088);
        keys.add(514);

        //keys.add(2615);
        //keys.add(3150);
        //keys.add(3750);
       // keys.add(898); keys.add(904); keys.add(1059); keys.add(1127); keys.add(1214);
       // keys.add(1219); keys.add(1220); keys.add(1243); keys.add(1250); keys.add(1252);
       // keys.add(1342);keys.add(1675);

        Collections.sort(keys);
        System.out.println("Total number of scans: " + keys.size());
        for (int scanId : keys) {
            Scan scan = scans.get(scanId);
            System.out.println("Analyzing scan " + scanId);
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            Collections.sort(peaks);
            GraphUtil.generateGapEdges(conf, peaks, 1);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            final Map<Integer, Integer> score = new HashMap<Integer, Integer>();
            System.out.println(scanId);
            List<String> tags = Tag4Finder.getTags(componentsFromGraph);
            for (String tag : tags) {
                System.out.println(tag);
            }
            Tag4Finder.fillScoresSeparator(componentsFromGraph, score, 300);
            if (!score.isEmpty()) {
                List<Integer> proteinIds = new ArrayList<Integer>();
                proteinIds.addAll(score.keySet());
                Collections.sort(proteinIds, new Comparator<Integer>() {
                    public int compare(Integer p1, Integer p2) {
                        return score.get(p2) - score.get(p1);
                    }
                });
                int[] pIds = new int[Math.min(20, proteinIds.size())];
                for (int i = 0; i < pIds.length; i++) {
                    pIds[i] = proteinIds.get(i);
                }       
                pIds = new int[]{3299};
                System.out.println("pIds.length = " + pIds.length);
                for (int i = 0; i < 3; i++) {
                    double[] bestPrsm = EValueAdapter.getBestPrsm(scan, pIds);

                    System.out.println(((int)bestPrsm[1]) +  " " + bestPrsm[0]);
                }
                //double oldEvalue = EValueAdapterOld.getBestEValue(scan, 3296);


                //System.out.println("old " + oldEvalue);
                //bestPrsm = EValueAdapter.getBestPrsm(scan, pIds);

                //System.out.println(proteinIds.get((int)bestPrsm[1]) +  " " + bestPrsm[0]);
                if (true) {
                    //break;
                }
                for (int i = 0; i < proteinIds.size(); i++) {
                    if (i == 20) {
                        break;
                    }
                    int proteinId = proteinIds.get(i);
                    double evalue = EValueAdapter.getBestPrsm(scan, proteinId)[0];
                    //oldEvalue = EValueAdapterOld.getBestEValue(scan, proteinId);
                    System.out.println("keyout " + proteinId + " " + evalue +  " " + score.get(proteinId) + " " + proteins.get(proteinId).getName());
                }
                System.out.println("Unchecked proteins:");
                for (int i = 20; i < proteinIds.size(); i++) {
                    int proteinId = proteinIds.get(i);
                    System.out.println((i + 1) + " " + proteinId + " " + score.get(proteinId));
                }

            }

            System.out.println("Finished scan analyzes");

        }


    }
}
