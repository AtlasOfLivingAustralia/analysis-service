/**
 * ************************************************************************
 * Copyright (C) 2010 Atlas of Living Australia All Rights Reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
 * the specific language governing rights and limitations under the License.
 * *************************************************************************
 */
package au.org.ala.spatial.util;

import au.org.ala.layers.intersect.Grid;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

/**
 * Class for region cutting test data grids
 *
 * @author adam
 */
public class GridCutterAloc {

    public static ArrayList<Object> loadCutGridsForAloc(File[] files, String extentsFilename, int pieces, AnalysisJob job) {
        ArrayList<Object> data = new ArrayList<Object>();

        if (job != null) {
            job.setProgress(0);
        }

        //determine outer bounds of layers
        double xmin = Double.MAX_VALUE;
        double ymin = Double.MAX_VALUE;
        double xmax = Double.MAX_VALUE * -1;
        double ymax = Double.MAX_VALUE * -1;
        double xres = 0.01;
        double yres = 0.01;
        for (File f : files) {
            String gridFilename = f.getPath().substring(0, f.getPath().length() - 4);
            Grid g = new Grid(gridFilename);
            xres = g.xres;
            yres = g.xres;
            if (xmin > g.xmin) {
                xmin = g.xmin;
            }
            if (xmax < g.xmax) {
                xmax = g.xmax;
            }
            if (ymin > g.ymin) {
                ymin = g.ymin;
            }
            if (ymax < g.ymax) {
                ymax = g.ymax;
            }
        }

        if (files.length < 2) {
            if (job != null) {
                job.setCurrentState(AnalysisJob.FAILED);
                job.log("Fewer than two layers with postive range.");
            } else {
                SpatialLogger.log("Fewer than two layers with postive range.");
            }
            return null;
        }


        //determine range and width's
        double xrange = xmax - xmin;
        double yrange = ymax - ymin;
        int width = (int) Math.ceil(xrange / xres);
        int height = (int) Math.ceil(yrange / yres);

        //write extents into a file now
        if (extentsFilename != null) {
            try {
                FileWriter fw = new FileWriter(extentsFilename);
                fw.append(String.valueOf(width)).append("\n");
                fw.append(String.valueOf(height)).append("\n");
                fw.append(String.valueOf(xmin)).append("\n");
                fw.append(String.valueOf(ymin)).append("\n");
                fw.append(String.valueOf(xmax)).append("\n");
                fw.append(String.valueOf(ymax));
                fw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        if (job != null) {
            job.setProgress(0.1, "exported extents");
        }

        //make cells list for outer bounds
        int th = height;
        int tw = width;
        int tp = 0;
        int[][] cells = new int[tw * th][2];
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                cells[tp][0] = j;
                cells[tp][1] = i;
                tp++;
            }
        }

        if (job != null) {
            job.setProgress(0.2, "determined target cells");
        }

        if (job != null) {
            job.log("Cut cells count: " + cells.length);
        } else {
            System.out.println("Cut cells count: " + cells.length);
        }

        //transform cells numbers to long/lat numbers
        double[][] points = new double[cells.length][2];
        for (int i = 0; i < cells.length; i++) {
            points[i][0] = xmin + cells[i][0] * xres;
            points[i][1] = ymin + cells[i][1] * yres;
        }

        //initialize data structure to hold everything
        // each data piece: row1[col1, col2, ...] row2[col1, col2, ...] row3...
        int remainingLength = cells.length;
        int step = (int) Math.floor(remainingLength / (double) pieces);
        for (int i = 0; i < pieces; i++) {
            if (i == pieces - 1) {
                data.add(new float[remainingLength * files.length]);
            } else {
                data.add(new float[step * files.length]);
                remainingLength -= step;
            }
        }

        //iterate for layers
        double[] layerExtents = new double[files.length * 2];
        for (int j = 0; j < files.length; j++) {
            String gridFilename = files[j].getPath().substring(0, files[j].getPath().length() - 4);
            Grid g = new Grid(gridFilename);
            float[] v = g.getValues2(points);

            //row range standardization
            float minv = Float.MAX_VALUE;
            float maxv = Float.MAX_VALUE * -1;
            for (int i = 0; i < v.length; i++) {
                if (v[i] < minv) {
                    minv = v[i];
                }
                if (v[i] > maxv) {
                    maxv = v[i];
                }
            }
            float range = maxv - minv;
            if (range > 0) {
                for (int i = 0; i < v.length; i++) {
                    v[i] = (v[i] - minv) / range;
                }
            } else {
                for (int i = 0; i < v.length; i++) {
                    v[i] = 0;
                }
            }
            layerExtents[j * 2] = minv;
            layerExtents[j * 2 + 1] = maxv;

            //iterate for pieces
            for (int i = 0; i < pieces; i++) {
                float[] d = (float[]) data.get(i);
                for (int k = j, n = i * step; k < d.length; k += files.length, n++) {
                    d[k] = v[n];
                }
            }

            if (job != null) {
                job.setProgress(0.2 + j / (double) files.length * 7 / 10.0, "opened grid: " + files[j].getName());
            }
        }

        if (job != null) {
            job.log("finished opening grids");
        }

        //remove null rows from data and cells
        int newCellPos = 0;
        int currentCellPos = 0;
        for (int i = 0; i < pieces; i++) {
            float[] d = (float[]) data.get(i);
            int newPos = 0;
            for (int k = 0; k < d.length; k += files.length) {
                int nMissing = 0;
                for (int j = 0; j < files.length; j++) {
                    if (Float.isNaN(d[k + j])) {
                        nMissing++;
                    }
                }
                //if (nMissing < files.length) {
                if (nMissing == 0) {
                    if (newPos < k) {
                        for (int j = 0; j < files.length; j++) {
                            d[newPos + j] = d[k + j];
                        }
                    }
                    newPos += files.length;
                    if (newCellPos < currentCellPos) {
                        cells[newCellPos][0] = cells[currentCellPos][0];
                        cells[newCellPos][1] = cells[currentCellPos][1];
                    }
                    newCellPos++;
                }
                currentCellPos++;
            }
            if (newPos < d.length) {
                d = java.util.Arrays.copyOf(d, newPos);
                data.set(i, d);
            }
        }

        //remove zero length data pieces
        for (int i = pieces - 1; i >= 0; i--) {
            float[] d = (float[]) data.get(i);
            if (d.length == 0) {
                data.remove(i);
            }
        }

        //add cells reference to output
        data.add(cells);

        //add extents to output
        double[] extents = new double[6 + layerExtents.length];
        extents[0] = width;
        extents[1] = height;
        extents[2] = xmin;
        extents[3] = ymin;
        extents[4] = xmax;
        extents[5] = ymax;
        for (int i = 0; i < layerExtents.length; i++) {
            extents[6 + i] = layerExtents[i];
        }
        data.add(extents);

        if (job != null) {
            job.setProgress(1, "cleaned data");
        }

        return data;
    }
}
