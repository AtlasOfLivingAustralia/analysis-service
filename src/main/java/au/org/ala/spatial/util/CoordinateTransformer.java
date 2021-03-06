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

import java.io.*;

/**
 * @author ajay
 */
public class CoordinateTransformer {

    public static String transformToGoogleMercator(String imgfilepath) {
        Runtime runtime = Runtime.getRuntime();
        try {

            String gdal_path = AlaspatialProperties.getGdalDir();
            //if (!gdal_path.endsWith("/")) gdal_path += "/";

            System.out.println("Got gdal_path: " + gdal_path);


            String base_command = gdal_path + "gdalwarp -s_srs EPSG:4326 -t_srs EPSG:900913 ";

            File oImg = new File(imgfilepath);


            String filenamepart = oImg.getName().substring(0, oImg.getName().lastIndexOf("."));
            String command = base_command + imgfilepath + " " + oImg.getParent() + File.separator + "t_" + filenamepart + ".tif";

            System.out.println("Exec'ing " + command);
            Process proc = runtime.exec(command);

            System.out.println("Setting up output stream readers");
            InputStreamReader isr = new InputStreamReader(proc.getInputStream());
            InputStreamReader eisr = new InputStreamReader(proc.getErrorStream());
            BufferedReader br = new BufferedReader(isr);
            BufferedReader ebr = new BufferedReader(eisr);
            String line;

            System.out.printf("Output of running %s is:", command);

            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }

            while ((line = ebr.readLine()) != null) {
                System.out.println(line);
            }

            int exitVal = proc.waitFor();

            // any error???
            //return exitVal;

            System.out.println(exitVal);

            // success == 0
            // if successful, then rename to the original filename
            if (exitVal == 0) {
                //oImg.delete();
                //File tImg = new File (oImg.getAbsolutePath() + "/t_" + oImg.getName());
                //tImg.renameTo(oImg);

                // now let's convert the image
                int cExitVal = convertImageType(oImg.getParent() + File.separator + "t_" + filenamepart + ".tif", oImg.getParent() + File.separator + "t_" + oImg.getName());

                System.out.println("Convert Image Type: " + cExitVal);

                return oImg.getParent() + File.separator + "t_" + oImg.getName();
            }
        } catch (Exception e) {
            System.out.println("OOOOPPPSSS: " + e.toString());
            System.out.println("{success: false , responseText: 'Error occurred' + " + e.toString() + "}");
            e.printStackTrace(System.out);
        }

        return null;
    }

    public static String transformAscToGeotif(String ascFile) {
        Runtime runtime = Runtime.getRuntime();
        try {

            String gdal_path = AlaspatialProperties.getGdalDir();

            System.out.println("Got gdal_path: " + gdal_path);

            String base_command = gdal_path + "gdal_translate -of GTiff -ot Float32 ";

            File oImg = new File(ascFile);

            String filenamepart = oImg.getName().substring(0, oImg.getName().lastIndexOf("."));
            String command = base_command + ascFile + " " + oImg.getParent() + File.separator + filenamepart + ".tif";

            System.out.println("Exec'ing " + command);
            Process proc = runtime.exec(command);

            System.out.println("Setting up output stream readers");
            InputStreamReader isr = new InputStreamReader(proc.getInputStream());
            InputStreamReader eisr = new InputStreamReader(proc.getErrorStream());
            BufferedReader br = new BufferedReader(isr);
            BufferedReader ebr = new BufferedReader(eisr);
            String line;

            System.out.printf("Output of running %s is:", command);

            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }

            while ((line = ebr.readLine()) != null) {
                System.out.println(line);
            }

            int exitVal = proc.waitFor();

            System.out.println(exitVal);

            if (exitVal == 0) {
                return oImg.getParent() + File.separator + oImg.getName() + ".tif";
            }
        } catch (Exception e) {
            System.out.println("OOOOPPPSSS: " + e.toString());
            System.out.println("{success: false , responseText: 'Error occurred' + " + e.toString() + "}");
            e.printStackTrace(System.out);
        }

        return null;
    }

    private static int convertImageType(String srcImgPath, String destImgPath) {
        try {

            Runtime runtime = Runtime.getRuntime();

            String command = AlaspatialProperties.getImageMagickDir() + " " + srcImgPath + " " + destImgPath;
            System.out.println("Exec'ing " + command);
            Process proc = runtime.exec(command);

            System.out.println("Setting up output stream readers");
            InputStreamReader isr = new InputStreamReader(proc.getInputStream());
            InputStreamReader eisr = new InputStreamReader(proc.getErrorStream());
            BufferedReader br = new BufferedReader(isr);
            BufferedReader ebr = new BufferedReader(eisr);
            String line;

            System.out.printf("Output of running %s is:", command);

            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }

            while ((line = ebr.readLine()) != null) {
                System.out.println(line);
            }

            int exitVal = proc.waitFor();

            // any error???
            //return exitVal;

            System.out.println(exitVal);

            return exitVal;

        } catch (Exception e) {
            System.out.println("Unable to convert image: ");
            e.printStackTrace(System.out);
        }

        return -1;
    }

    public static void generateWorldFiles(String outputpath, String baseFilename, String xRes, String yRes, String xMin, String yMin) {
        try {

            if (!outputpath.endsWith(File.separator)) {
                outputpath += File.separator;
            }

            System.out.println("Generating world files for " + baseFilename + " under " + outputpath);

            StringBuffer sbWorldFile = new StringBuffer();
            //  pixel X size
            // sbWorldFile.append(xRes).append("\n");
            sbWorldFile.append(xRes).append("\n");
            // rotation about the Y axis (usually 0.0)
            sbWorldFile.append("0").append("\n");
            // rotation about the X axis (usually 0.0)
            sbWorldFile.append("0").append("\n");
            // negative pixel Y size
            // sbWorldFile.append(yRes).append("\n");
            sbWorldFile.append(yRes).append("\n");
            // X coordinate of upper left pixel center
            // sbWorldFile.append(xMin).append("\n");
            sbWorldFile.append(xMin).append("\n");
            // Y coordinate of upper left pixel center
            // sbWorldFile.append(yMin).append("\n");
            sbWorldFile.append(yMin).append("\n");

            PrintWriter pgwout = new PrintWriter(new BufferedWriter(new FileWriter(outputpath + baseFilename + ".pgw")));
            pgwout.write(sbWorldFile.toString());
            pgwout.close();

            generate4326prj(outputpath, baseFilename);

        } catch (Exception e) {
            e.printStackTrace(System.out);
        }

    }

    public static void generate4326prj(String outputpath, String baseFilename) {
        try {

            if (!outputpath.endsWith(File.separator)) {
                outputpath += File.separator;
            }

            System.out.println("Generating prj file for " + baseFilename + " under " + outputpath);

            StringBuffer sbProjection = new StringBuffer();
            sbProjection.append("GEOGCS[\"WGS 84\", ").append("\n");
            sbProjection.append("    DATUM[\"WGS_1984\", ").append("\n");
            sbProjection.append("        SPHEROID[\"WGS 84\",6378137,298.257223563, ").append("\n");
            sbProjection.append("            AUTHORITY[\"EPSG\",\"7030\"]], ").append("\n");
            sbProjection.append("        AUTHORITY[\"EPSG\",\"6326\"]], ").append("\n");
            sbProjection.append("    PRIMEM[\"Greenwich\",0, ").append("\n");
            sbProjection.append("        AUTHORITY[\"EPSG\",\"8901\"]], ").append("\n");
            sbProjection.append("    UNIT[\"degree\",0.01745329251994328, ").append("\n");
            sbProjection.append("        AUTHORITY[\"EPSG\",\"9122\"]], ").append("\n");
            sbProjection.append("    AUTHORITY[\"EPSG\",\"4326\"]] ").append("\n");

            PrintWriter prjout = new PrintWriter(new BufferedWriter(new FileWriter(outputpath + baseFilename + ".prj")));
            prjout.write(sbProjection.toString());
            prjout.close();
        } catch (Exception e) {
            e.printStackTrace(System.out);
        }
    }

    public static String transformAscToGeotiff(String src) {
        Runtime runtime = Runtime.getRuntime();
        try {

            String gdal_path = AlaspatialProperties.getGdalDir();

            System.out.println("Got gdal_path: " + gdal_path);

            String base_command = gdal_path + "gdal_translate -of GTiff ";

            String command = base_command + src + " " + src + ".tif";

            System.out.println("Exec'ing " + command);
            Process proc = runtime.exec(command);

            System.out.println("Setting up output stream readers");
            InputStreamReader isr = new InputStreamReader(proc.getInputStream());
            InputStreamReader eisr = new InputStreamReader(proc.getErrorStream());
            BufferedReader br = new BufferedReader(isr);
            BufferedReader ebr = new BufferedReader(eisr);
            String line;

            System.out.printf("Output of running %s is:", command);

            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }

            while ((line = ebr.readLine()) != null) {
                System.out.println(line);
            }

            int exitVal = proc.waitFor();

            // any error???
            //return exitVal;

            System.out.println(exitVal);

            // success == 0
            // if successful, then rename to the original filename
            if (exitVal == 0) {
                return src + ".tif";
            }
        } catch (Exception e) {
            System.out.println("OOOOPPPSSS: " + e.toString());
            System.out.println("{success: false , responseText: 'Error occurred' + " + e.toString() + "}");
            e.printStackTrace(System.out);
        }

        return null;
    }

    public static void diva2asc(String diva, String asc) {
        Grid g = new Grid(diva);
        float[] grid_data = g.getGrid();

        //export ASCGRID
        BufferedWriter fw = null;
        try {
            fw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(asc + ".asc"), "US-ASCII"));
            fw.append("ncols ").append(String.valueOf(g.ncols)).append("\n");
            fw.append("nrows ").append(String.valueOf(g.nrows)).append("\n");
            fw.append("xllcorner ").append(String.valueOf(g.xmin)).append("\n");
            fw.append("yllcorner ").append(String.valueOf(g.ymin)).append("\n");
            fw.append("cellsize ").append(String.valueOf(g.xres)).append("\n");

            fw.append("NODATA_value ").append(String.valueOf(-1));

            for (int i = 0; i < g.nrows; i++) {
                fw.append("\n");
                for (int j = 0; j < g.ncols; j++) {
                    if (j > 0) {
                        fw.append(" ");
                    }
                    if (Double.isNaN(grid_data[i * g.ncols + j])) {
                        fw.append("-1");
                    } else {
                        fw.append(String.valueOf(grid_data[i * g.ncols + j]));
                    }
                }
            }
            fw.append("\n");
        } catch (Exception e) {
            e.printStackTrace(System.out);
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (Exception e) {
                    e.printStackTrace(System.out);
                }
            }
        }
    }
}
