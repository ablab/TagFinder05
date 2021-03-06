package ru.spbau.bioinf.tagfinder;

import java.text.NumberFormat;
import java.util.*;

public class FourSearch {
    private static NumberFormat df = NumberFormat.getInstance();
    private static List<Protein> proteins;
    private static Configuration conf;

    static {
        df.setMaximumFractionDigits(2);
    }

    public static void main(String[] args) throws Exception {
        long start = System.currentTimeMillis();                
        conf = new Configuration(args);
        EValueAdapter.init(conf, "TARGET");
        proteins = conf.getProteins();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        for (Scan scan : scans.values()) {
            if (scan.getPeaks().size() >= 10 && scan.getPrecursorMass() > 2500) {
                keys.add(scan.getId());
            }
        }
        Collections.sort(keys);
        Set<Integer> forProcess = new HashSet<Integer>();
        forProcess.addAll(conf.getBadMSAlignResults().keySet());
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

        //keys.clear(); keys.add(1250);

        System.out.println("Total number of scans: " + keys.size());
        for (int scanId : keys) {
            if (msAlignResults.keySet().contains(scanId)) {
                //continue;
            }
            if (!forProcess.contains(scanId)) {
                //continue;
            }
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            Collections.sort(peaks);
            //GraphUtil.generateEdges(conf, peaks);
            GraphUtil.generateGapEdges(conf, peaks, 1);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            final Map<Integer, Integer> score = new HashMap<Integer, Integer>();
            System.out.print(scanId) ;
            for (List<Peak> component : componentsFromGraph) {
                Peak[] tagPeaks = GraphUtil.findBestTag(component);
                /*
                int len = tagPeaks.length;
                if (len > 4) {
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
                if (tagPeaks.length < 4) {
                    continue;
                }
                for (int j = 0; j < tagPeaks.length  - 1; j++) {
                    superTag += Acid.getAcid(tagPeaks[j + 1].getValue() - tagPeaks[j].getValue()).name();
                }
                //superTag = "AVSFMV";
                //System.out.print(superTag + " ");

                for (int i = 0; i < superTag.length() - 3; i++) {
                    String tag = superTag.substring(i, i + 4);
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
                double bestEvalue = 10E9;
                int bestProteinId = -1;
                for (int i = 0; i < proteinIds.size(); i++) {
                    if (i == 20) {
                        break;
                    }
                    int proteinId = proteinIds.get(i);
                    if (score.get(proteinId) > 0) {
                        //System.out.println(proteinId + " " + score.get(proteinId) + " " + proteins.get(proteinId).getName());
                    }

                    double evalue = EValueAdapter.getBestPrsm(scan, proteinId)[0];
                    if (evalue < 1) {
                        if (evalue < Configuration.EVALUE_LIMIT) {
                            //System.out.print("!");
                        }
                        if (evalue < bestEvalue) {
                            bestProteinId = proteinId;
                            bestEvalue = evalue;
                        }
                        //System.out.println("Success: " + evalue + " " +  proteinId + " " + score.get(proteinId) + " " + proteins.get(proteinId).getName());
                    }
                }
                if (bestProteinId >= 0) {
                    System.out.print(" " + bestProteinId + " " + bestEvalue +  " " +  score.get(bestProteinId) + " " + proteins.get(bestProteinId).getName());
                }
            }
            System.out.println();

        }

        System.out.println("time: " + (System.currentTimeMillis() - start));
    }
}
