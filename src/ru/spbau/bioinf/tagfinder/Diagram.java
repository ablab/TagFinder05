package ru.spbau.bioinf.tagfinder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import ru.spbau.bioinf.tagfinder.util.ReaderUtil;

public class Diagram {
    public static void main(String[] args) throws Exception {
        Set<Integer> matchedNew = new HashSet<Integer>();
        File dir = new File("");
        dir = new File(dir.getCanonicalPath());
        File[] txtFiles = dir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith("out20.txt");
            }
        });
        int lines = 0;
        for (File txt : txtFiles) {
            BufferedReader in = ReaderUtil.createInputReader(txt);
            String s = in.readLine();
            if (!s.startsWith("#Total ")) {
                continue;
            }
            System.out.println("Processing file " + txt.getName());

            while ((s = in.readLine())!= null) {
                if (!s.startsWith("#")) {
                    lines++;
                    String[] data = s.split(" ");
                    if (data.length > 1) {
                        if (Double.parseDouble(data[2]) < 0.0024) {
                            matchedNew.add(Integer.parseInt(data[0]));
                        }
                    }
                }
            }
        }
        System.out.println(lines + " results processed.");
        BufferedReader readerMsalign = ReaderUtil.createInputReader(new File("msalign.txt"));
        Set<Integer> msAlignOld = new HashSet<Integer>();
        String s;
        while ((s = readerMsalign.readLine())!= null) {
            msAlignOld.add(Integer.parseInt(s));
        }

        Integer[] mn = matchedNew.toArray(new Integer[]{});
        Arrays.sort(mn);
        for (Integer integer : mn) {
            System.out.println(integer);
        }



        Configuration conf = new Configuration(args);
        Map<Integer, Scan> scans = conf.getScans();
        List<Integer> keys = new ArrayList<Integer>();
        for (Scan scan : scans.values()) {
            if (scan.getPeaks().size() >= 10 && scan.getPrecursorMass() > 2500) {
                keys.add(scan.getId());
            }
        }
        Collections.sort(keys);


        int msAlign = 0;
        int tag = 0;
        int both = 0;
        int total = 0;
        int tagOnly = 0;
        int msOnly = 0;
        for (int scanId : keys) {
            total++;
            boolean isNew = matchedNew.contains(scanId);
            boolean isOld = msAlignOld.contains(scanId);
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
                System.out.println("New " + scanId);
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
