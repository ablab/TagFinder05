package ru.spbau.bioinf.tagfinder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FindSubtagInDatabase {

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        EValueAdapter.init(conf, "TARGET");
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        keys.add(666); keys.add(906); keys.add(1038); keys.add(1082); keys.add(1088); keys.add(1186); keys.add(1251); keys.add(1252);
        keys.add(1256); keys.add(1322); keys.add(1608); keys.add(1715); keys.add(1718); keys.add(2180); keys.add(2480);
        Collections.sort(keys);
        List<Protein> proteins = conf.getProteins();
        Tag4Finder.proteins = proteins;
        for (int scanId : keys) {
            System.out.println("Analyzing scan #" + scanId);
            Scan scan = scans.get(scanId);
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            GraphUtil.generateEdges(conf, peaks);
            Set<String> tags = GraphUtil.generateTags(conf, peaks);
            int maxLen = 0;
            Set<String> process = new HashSet<String>();
            for (String tag : tags) {
                if (tag.length() >= maxLen) {
                    for (Protein protein : proteins) {
                        if (protein.getSimplifiedAcids().contains(tag)) {
                            process.add(tag);
                            maxLen = tag.length();
                        }
                    }
                }
            }

            for (Iterator<String> iterator = process.iterator(); iterator.hasNext(); ) {
                String next = iterator.next();
                if (next.length() < maxLen) {
                    iterator.remove();
                }
            }
            System.out.println("maxLen = " + maxLen);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            final Map<Integer, Integer> score = new HashMap<Integer, Integer>();
            Tag4Finder.fillScoresSeparator(componentsFromGraph, score, 300);

            for (String tag : process) {
                System.out.println(tag);
                for (Protein protein : proteins) {
                    if (protein.getSimplifiedAcids().contains(tag)) {
                        int proteinId = protein.getProteinId();
                        System.out.println(proteinId + " " + score.get(proteinId) + " " + EValueAdapter.getBestPrsm(scan, proteinId)[0]);
                    }
                }
            }
        }
    }
}
