package ru.spbau.bioinf.tagfinder;


import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.plugins.jpeg.JPEGImageWriteParam;
import javax.imageio.stream.ImageOutputStream;

public class ImageGenerator {

    public static DecimalFormat df = (DecimalFormat) NumberFormat.getNumberInstance();

    static {
        df.setMaximumFractionDigits(2);
        DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        df.setDecimalFormatSymbols(dfs);
        df.setGroupingSize(10000);
    }

    private static Color[] colors = new Color[] {Color.RED, Color.GREEN, Color.MAGENTA, Color.BLUE, Color.ORANGE};
    
    public static void main(String[] args) throws Exception {    
        Configuration conf = new Configuration(args);
        Map<Integer, Scan> scans = conf.getScans();
        File imagesDir = new File("images");
        imagesDir.mkdirs();
        for (Scan scan : scans.values()) {
            List<Peak> peaks = scan.createSpectrumWithYPeaks(PrecursorMassShiftFinder.getPrecursorMassShiftForMoreEdges(conf, scan));
            Collections.sort(peaks);
            GraphUtil.generateGapEdges(conf, peaks, 1);
            List<List<Peak>> componentsFromGraph = GraphUtil.getComponentsFromGraph(peaks);
            List<List<Peak[]>> a = new ArrayList<List<Peak[]>>();
            int totalLines = 0;
            for (List<Peak> component : componentsFromGraph) {
                if (component.size() >= 4) {
                    List<Peak[]> lines = new ArrayList<Peak[]>();
                    do {
                        Peak[] bestTag = GraphUtil.findBestTag(component);
                        if (bestTag.length < 2) {
                            break;
                        }
                        lines.add(bestTag);
                        for (int i = 0; i < bestTag.length - 1; i++) {
                            bestTag[i].removeNext(bestTag[i+1]);
                        }
                    } while (true);
                    totalLines += lines.size();
                    a.add(lines);
                }                
            }

            totalLines += a.size();

            if (totalLines == 0) {
                continue;
            }
            int lineHeight = 25;
            int width = 1000;
            int height = totalLines * lineHeight;
            BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR );
            Graphics g = bi.getGraphics();
            g.setColor(Color.WHITE);
            g.fillRect(0, 0, width, height);
            g.setColor(Color.BLACK);

            int curLine = 0;
            for (int componentId = 1; componentId <= a.size(); componentId++) {
                List<Peak[]> lines = a.get(componentId - 1);
                double min = 10E9;
                Set<Peak> used = new HashSet<Peak>();
                Set<Peak> forColor = new HashSet<Peak>();
                
                for (Peak[] line : lines) {
                    for (Peak peak : line) {
                        min = Math.min(peak.getValue(), min);
                        if (used.contains(peak)) {
                            forColor.add(peak);
                        } else {
                            used.add(peak);
                        }
                    }
                }

                List<Peak> duplicates = new ArrayList<Peak>();
                duplicates.addAll(forColor);
                Collections.sort(duplicates);
                Map<Peak, Integer> colorMap = new HashMap<Peak, Integer>();
                for (int i = 0; i < duplicates.size(); i++) {
                    colorMap.put(duplicates.get(i), i);                    
                }

                g.setColor(Color.BLACK);
                g.drawString("Component " + componentId + " starting from " + df.format(min), 5, (curLine + 1)* lineHeight - 3);
                curLine++;
                
                for (Peak[] line : lines) {
                    for (int i = 0; i < line.length; i++) {
                        Peak peak = line[i];
                        int x = (int)Math.round((peak.getValue() - min) /3 + 5);
                        int y1 = curLine * lineHeight + 3;
                        int y2 = (curLine + 1) * lineHeight - 3;
                        Color color = Color.BLACK;
                        if (colorMap.containsKey(peak)) {
                            color = colors[colorMap.get(peak) % colors.length];
                        }
                        g.setColor(color);
                        g.drawLine(x, y1, x, y2);
                        int dir = peak.getPeakType() == PeakType.B ? -1 : 1;
                        g.drawLine(x + dir * 3, y1, x, y1);
                        g.drawLine(x + dir * 3, y2, x, y2);
                        if (i < line.length - 1) {
                            double m1 = peak.getValue();
                            double m2 = line[i + 1].getValue();
                            Acid acid = Acid.getAcid(m2 - m1);
                            int x2 = (int)Math.round((m2 - min) /3 + 5);
                            int c = (x + x2)/2;
                            g.setColor(Color.BLACK);
                            g.drawString(acid.name(), c - 5, (curLine + 1)* lineHeight - 7);
                        }

                    }
                    curLine++;
                }
            }

            ImageWriter writer = ImageIO.getImageWritersByFormatName("jpg").next();
            ImageOutputStream ios = ImageIO.createImageOutputStream( new File(imagesDir, "spectrum" + scan.getId() + ".jpg"));
            writer.setOutput(ios);
            ImageWriteParam param = new JPEGImageWriteParam( java.util.Locale.getDefault() );
            param.setCompressionMode(ImageWriteParam.MODE_EXPLICIT) ;
            param.setCompressionQuality(0.98f);

            writer.write(null, new IIOImage(bi, null, null ), param);
        }
    }
}
