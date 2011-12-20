package ru.spbau.bioinf.tagfinder;

import edu.ucsd.msalign.spec.id.EValueAdapter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FourSearchPositive {
    private static NumberFormat df = NumberFormat.getInstance();
    private static List<Protein> proteins;
    private static Configuration conf;

    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        conf = new Configuration(args);
        EValueAdapter.init(conf);
        proteins = conf.getProteins();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(scans.keySet());
        Collections.sort(keys);
        Set<Integer> forProcess = new HashSet<Integer>();
        forProcess.addAll(msAlignResults.keySet());
        /*

        forProcess.add(1088);
        forProcess.add(1243);
        forProcess.add(1247);
        forProcess.add(1250);
        forProcess.add(1251);
        forProcess.add(1252);
        forProcess.add(1675);
        forProcess.add(1718);
        forProcess.add(1726);
        forProcess.add(2239);
        forProcess.add(3695);

        forProcess.clear();forProcess.add(1247);
        */

        int good = 0;
        int bad = 0;

        int TAG_SIZE = 3;

        for (int scanId : keys) {
            if (!forProcess.contains(scanId)) {
                continue;
            }
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            Collections.sort(peaks);
            //GraphUtil.generateEdges(conf, peaks);
            GraphUtil.generateGapEdges(conf, peaks, 1);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            final Map<Integer, Integer> score = new HashMap<Integer, Integer>();
            System.out.println("Scan " + scanId) ;
            for (List<Peak> component : componentsFromGraph) {
                Peak[] tagPeaks = GraphUtil.findBestTag(component);/*
                int len = tagPeaks.length;
                if (len > 5) {
                    Collections.sort(component, new Comparator<Peak>() {
                        public int compare(Peak o1, Peak o2) {
                            double intensity1 = o1.getIntensity();
                            double intensity2 = o2.getIntensity();
                            if (intensity1 == intensity2) {
                                return 0;
                            }
                            return intensity2 - intensity1 > 0 ? 1 : -1;
                        }
                    });
                    while (true) {
                        component.remove(component.size() - 1);
                        for (Peak peak : component) {
                            peak.clearEdges();
                        }
                        GraphUtil.generateEdges(conf, peaks);
                        Peak[] tagPeaks2 = GraphUtil.findBestTag(component);
                        if (tagPeaks2.length == tagPeaks.length) {
                            tagPeaks = tagPeaks2;
                        } else {
                            break;
                        }
                    }
                }
                */
                String superTag = "";
                if (tagPeaks.length < TAG_SIZE  + 1) {
                    continue;
                }
                for (int j = 0; j < tagPeaks.length  - 1; j++) {
                    superTag += Acid.getAcid(tagPeaks[j + 1].getValue() - tagPeaks[j].getValue()).name();
                }
                //superTag = "AVSFMV";
                //System.out.print(superTag + " ");

                for (int i = 0; i < superTag.length() - TAG_SIZE + 1; i++) {
                    String tag = superTag.substring(i, i + TAG_SIZE);
                    for (Protein protein : proteins) {
                        if (protein.getSimplifiedAcids().contains(tag)) {
                            int proteinId = protein.getProteinId();
                            if (!score.containsKey(proteinId)) {
                                score.put(proteinId, 0);
                            }
                            score.put(proteinId, score.get(proteinId) + 1);
                        }
                    }
                }
            }
            if (!score.isEmpty()) {
                List<Integer> proteinIds = new ArrayList<Integer>();
                proteinIds.addAll(score.keySet());
                Collections.sort(proteinIds, new Comparator<Integer>() {
                    public int compare(Integer p1, Integer p2) {
                        return score.get(p2) - score.get(p1);
                    }
                });
                boolean hasGood = false;
                for (int i = 0; i < proteinIds.size(); i++) {
                    if (i == 20) {
                        break;
                    }
                    int proteinId = proteinIds.get(i);
                    if (score.get(proteinId) > 3) {
                        //System.out.println(proteinId + " " + score.get(proteinId) + " " + proteins.get(proteinId).getName());
                        if (msAlignResults.get(scanId) == proteinId) {
                            hasGood = true;
                        }
                    }
                    /*
                    double evalue = EValueAdapter.getBestEValue(scan, proteinId);
                    if (evalue < 100) {
                        System.out.println();
                        System.out.println(evalue);
                    }
                    */
                }
                if (hasGood) {
                    good++;
                } else {
                    bad++;
                }
            } else {
                bad++;
            }
            System.out.println("good " + good + " bad " + bad);
        }
    }
}
