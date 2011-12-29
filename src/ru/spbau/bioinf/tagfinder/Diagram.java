package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class Diagram {
    public static void main(String[] args) throws Exception {
        Set<Integer> matchedNew = new HashSet<Integer>();
        BufferedReader in = ReaderUtil.createInputReader(new File("tag4first.log"));
        String s;
        while ((s = in.readLine())!= null) {
            if (!s.startsWith("#")) {
                String[] data = s.split(" ");
                if (data.length > 1) {
                    if (Double.parseDouble(data[2]) < 0.0024) {
                        matchedNew.add(Integer.parseInt(data[0]));
                    }
                }
            }
        }
        Configuration conf = new Configuration(args);
        Map<Integer, Scan> scans = conf.getScans();
        Map<Integer, Integer> msAlignResults = conf.getMSAlignResults();
        Map<Integer, Double> evalues = conf.getEvalues();
        int msAlign = 0;
        int tag = 0;
        int both = 0;
        int total = 0;
        int tagOnly = 0;
        int msOnly = 0;
        List<Integer> keys = new ArrayList<Integer> ();
        keys.addAll(msAlignResults.keySet());
        Collections.sort(keys);
        for (int scanId : keys) {
            total++;
            Scan scan = scans.get(scanId);
            boolean isNew = matchedNew.contains(scanId);
            boolean isOld = evalues.get(scanId) < 0.0024;
            if (isOld) {
                //System.out.println(scanId);
            }
            if (isNew) {
                tag++;
            }
            if (isOld) {
                msAlign++;
            }
            if (isNew && isOld) {
                both++;
            }

            if (isNew && !isOld) {
                tagOnly++;
            }

            if (!isNew && isOld) {
                msOnly++;
            }
        }
        System.out.println("total = " + total);
        System.out.println("msAlign = " + msAlign);
        System.out.println("tag = " + tag);
        System.out.println("both = " + both);
        System.out.println("msOnly = " + msOnly);
        System.out.println("tagOnly = " + tagOnly);
    }
}
