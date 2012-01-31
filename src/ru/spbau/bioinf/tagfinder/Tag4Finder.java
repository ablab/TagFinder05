package ru.spbau.bioinf.tagfinder;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Tag4Finder {
    private static NumberFormat df = NumberFormat.getInstance();
    protected static List<Protein> proteins;
    private static Configuration conf;

    static {
        df.setMaximumFractionDigits(2);
    }

    static int TAG_SIZE = 4;
    
    public static void main(String[] args) throws Exception {
        String[] parameters = new String[]{};
        int start = -1;
        int candidatesCount = 20;
        int batchSize = candidatesCount;
        

        for (String arg : args) {
            int ind = arg.indexOf('=');
            if (ind > 0) {
                String parameter = arg.substring(0, ind);
                String value = arg.substring(ind + 1);
                if ("dataset".equals(parameter)) {
                    File file = new File(value);
                    if (file.isDirectory()) {
                        parameters = new String[]{value};
                    }
                }
                if ("start".equals(parameter)) {
                    try {
                        start = Integer.parseInt(arg);
                    } catch (NumberFormatException e) {
                    }
                }
                if ("batchSize".equals(parameter)) {
                    try {
                        batchSize = Integer.parseInt(value);
                    } catch (NumberFormatException e) {
                    }
                }
            }
        }
        System.out.println("candidatesCount = " + candidatesCount);
        System.out.println("batchSize = " + batchSize);
                
        conf = new Configuration(parameters);
        EValueAdapter.init(conf);
        proteins = conf.getProteins();
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        for (Scan scan : scans.values()) {
            if (scan.getPeaks().size() >= 10 && scan.getPrecursorMass() > 2500) {
                int scanId = scan.getId();
                if (scanId >= start) {
                    keys.add(scanId);
                }
            }
        }
        Collections.sort(keys);
        System.out.println("#Total number of scans for process: " + keys.size());
        System.out.println("#Started at " + new Date());
        for (int scanId : keys) {
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            Collections.sort(peaks);
            GraphUtil.generateGapEdges(conf, peaks, 1);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            final Map<Integer, Integer> score = new HashMap<Integer, Integer>();
            System.out.print(scanId + " ") ;
            fillScoresSeparator(componentsFromGraph, score);
            if (!score.isEmpty()) {
                List<Integer> proteinIds = new ArrayList<Integer>();
                proteinIds.addAll(score.keySet());
                Collections.sort(proteinIds, new Comparator<Integer>() {
                    public int compare(Integer p1, Integer p2) {
                        return score.get(p2) - score.get(p1);
                    }
                });
                double bestEvalue = 10E99;
                int bestProteinId = -1;
                Set<Integer> p = new HashSet<Integer>();
                
                for (int batchStart = 0; batchStart < candidatesCount; batchStart += batchSize) {
                    for (int i = batchStart; i < batchStart + batchSize; i++) {
                        if (i >= proteinIds.size()) {
                            break;
                        }
                        int proteinId = proteinIds.get(i);
                        p.add(proteinId);
                    }
                    int[] pIds = new int[p.size()];
                    int cur = 0;
                    for (Integer pId : p) {
                        pIds[cur] = pId;
                        cur++;
                    }
                    double[] best = EValueAdapter.getBestPrsm(scan, pIds);
                    if (best[1] > -1 && best[0] < bestEvalue) {
                        bestProteinId = pIds[(int) Math.round(best[1])];
                        bestEvalue = best[0];
                    }
                }

                if (bestProteinId >= 0) {
                    System.out.print(bestProteinId + " " + bestEvalue +  " " +  score.get(bestProteinId) + " " + proteins.get(bestProteinId).getName());
                }
            }
            System.out.println();
        }
        System.out.println("#Finished at " + new Date()) ;
    }


    public static void fillScoresSeparator(List<List<Peak>> componentsFromGraph, Map<Integer, Integer> scores) {
        List<String> tags = getTags(componentsFromGraph);
        int subSequenceLength = 300;        
        for (Protein protein : proteins) {
            int score = 0;
            String sequence = protein.getSimplifiedAcids();
            int end = sequence.length() - subSequenceLength + 1;
            if (end < 1) {
                end = 1;
            }
            for (int start = 0; start < end; start++) {
                String s = sequence.substring(start, Math.min(start + subSequenceLength, sequence.length()));
                int nextScore = 0;
                for (String tag : tags) {
                    for (int i = 0; i <= tag.length() - TAG_SIZE; i++) {
                        if (s.contains(tag.substring(i, i + TAG_SIZE))) {
                            nextScore++;
                        }
                    }
                }
                if (nextScore > score) {
                    score = nextScore;
                }
            }
            if (score > 0) {
                scores.put(protein.getProteinId(), score);
            }
        }
    }

    public static List<String> getTags(List<List<Peak>> componentsFromGraph) {
        List<String> tags = new ArrayList<String>();
        for (List<Peak> component : componentsFromGraph) {
            Peak[] tagPeaks = GraphUtil.findBestTag(component);
            String superTag = "";
            if (tagPeaks.length < TAG_SIZE  + 1) {
                continue;
            }
            for (int j = 0; j < tagPeaks.length  - 1; j++) {
                superTag += Acid.getAcid(tagPeaks[j + 1].getValue() - tagPeaks[j].getValue()).name();
            }
            tags.add(superTag);
        }
        return tags;
    }
}
