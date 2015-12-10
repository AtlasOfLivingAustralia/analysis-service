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

import au.org.ala.layers.client.Client;
import au.org.ala.layers.dao.LayerDAO;
import au.org.ala.layers.dto.Layer;
import au.org.ala.layers.grid.GridCutter;
import au.org.ala.layers.intersect.Grid;
import au.org.ala.layers.intersect.SimpleRegion;
import au.org.ala.layers.util.LayerFilter;
import au.org.ala.spatial.analysis.maxent.MaxentServiceImpl;
import au.org.ala.spatial.analysis.maxent.MaxentSettings;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.Hashtable;

/**
 * @author Adam
 */
public class AnalysisJobMaxent extends AnalysisJob {

    long[] stageTimes;
    String currentPath;
    String taxon;
    String area;
    String envlist;
    String txtTestPercentage;
    String chkJackknife;
    String chkResponseCurves;
    LayerFilter[] envelope;
    SimpleRegion region;
    int cells;
    int speciesCount;
    String resolution;
    String[] envnameslist;

    public AnalysisJobMaxent(String pid, String currentPath_, String taxon_, String envlist_, SimpleRegion region_, LayerFilter[] filter_, String txtTestPercentage_, String chkJackknife_, String chkResponseCurves_, String resolution_) {
        super(pid);
        currentPath = currentPath_;
        taxon = taxon_;
        region = region_;
        envelope = filter_;
        txtTestPercentage = txtTestPercentage_;
        chkJackknife = chkJackknife_;
        chkResponseCurves = chkResponseCurves_;
        envlist = envlist_;
        resolution = resolution_;
        envnameslist = envlist.split(":");
        System.out.println("ENVLIST*** " + envlist);


        //TODO: remove rough estimate
        if (region != null) {
            cells = (int) Math.ceil(region.getWidth() / Double.parseDouble(resolution)
                    * region.getHeight() / Double.parseDouble(resolution));
        } else {
            cells = 1000000; //or something
        }

        speciesCount = 10000;

        stageTimes = new long[4];
        setStage(0);
        setProgress(0);
    }

    static public void readReplace(String fname, String oldPattern, String replPattern) {
        String line;
        StringBuffer sb = new StringBuffer();
        try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
            while ((line = reader.readLine()) != null) {
                line = line.replaceAll(oldPattern, replPattern);
                sb.append(line + "\n");
            }
            reader.close();
            BufferedWriter out = new BufferedWriter(new FileWriter(fname));
            out.write(sb.toString());
            out.close();
        } catch (Throwable e) {
            e.printStackTrace(System.out);
        }
    }

    static public void readReplaceAfter(String fname, String start, String oldPattern, String replPattern) {
        String line;
        StringBuffer sb = new StringBuffer();
        try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
            int afterPos = -1;
            while ((line = reader.readLine()) != null) {
                if (afterPos < 0 && (afterPos = line.indexOf(start)) >= 0) {
                    line = line.substring(0, afterPos + start.length()) + line.substring(afterPos + start.length()).replaceAll(oldPattern, replPattern);
                } else if (afterPos > 0) {
                    line = line.replaceAll(oldPattern, replPattern);
                }
                sb.append(line + "\n");
            }
            reader.close();
            BufferedWriter out = new BufferedWriter(new FileWriter(fname));
            out.write(sb.toString());
            out.close();
        } catch (Throwable e) {
            e.printStackTrace(System.out);
        }
    }

    @Override
    public void run() {

        try {
            setCurrentState(RUNNING);
            setStage(0);

            // dump the species data to a file
            setProgress(0, "dumping species data");

            setProgress(0, "preparing input files and run parameters");

            //get species query data
            File speciesqueryfile = new File(currentPath + "output" + File.separator + "maxent" + File.separator + getName() + File.separator + "species_query.csv");
            BufferedReader br = new BufferedReader(new FileReader(speciesqueryfile));
            String speciesq = br.readLine();
            String bs = br.readLine();
            br.close();
            FileUtils.deleteQuietly(speciesqueryfile);
            OccurrenceData od = new OccurrenceData();
            String[] s = od.getSpeciesData(speciesq, bs, null, null);

            if (s[0] == null) {
                //error, no species
                System.out.println("Has error, sending maxent error message");
                setProgress(1, "failed: ");
                setCurrentState(FAILED);
                setMessage("No species selected.\nHint: Make sure your active area includes species occurrences");
            } else {

                writeFile(s[0], currentPath + "output" + File.separator + "maxent" + File.separator + getName() +
                        File.separator, "species_points.csv");
                if (s[1] != null) {
                    writeFile(s[1], currentPath + "output" + File.separator + "maxent" + File.separator + getName() +
                            File.separator, "Prediction_removedSpecies.txt");
                }

                //add layer display names to cutDataPath
                String names = "";
                String[] displayNames = new String[envnameslist.length];
                for (int i = 0; i < envnameslist.length; i++) {
                    String[] name_displayname = envnameslist[i].split("\\|");
                    if (name_displayname.length > 1) {
                        envnameslist[i] = name_displayname[0];
                        displayNames[i] = name_displayname[1] + " (" + envnameslist[i] + ")";
                        names += "\n" + envnameslist[i] + "=" + name_displayname[1];
                    } else {
                        envnameslist[i] = envnameslist[i];
                        displayNames[i] = envnameslist[i];
                        names += "\n" + envnameslist[i] + "=" + envnameslist[i];
                    }
                }

                String cutDataPath = GridCutter.cut2(envnameslist, resolution, region, envelope, null);

                FileUtils.writeStringToFile(new File(cutDataPath + File.separator + "additional_properties.txt"), names);

                MaxentSettings msets = new MaxentSettings();
                msets.setMaxentPath(AlaspatialProperties.getAnalysisMaxentCmd());
                msets.setRandomTestPercentage(Integer.parseInt(txtTestPercentage));
                msets.setEnvPath(cutDataPath);          //use (possibly) cut layers

                String ctxVarToggler = "";
                for (int l = 0; l < envnameslist.length; l++) {
                    if (Client.getLayerDao().getLayerByName(envnameslist[l]) != null
                            && Client.getLayerDao().getLayerByName(envnameslist[l]).getType().equals("contextual")) {
                        ctxVarToggler += envnameslist[l] + " ";
                    } else if (envnameslist[l].startsWith("aloc_")) {
                        ctxVarToggler += envnameslist[l] + " ";
                    }
                }
                msets.setEnvVarToggler(ctxVarToggler);
                msets.setEnvList(Arrays.asList(envnameslist.clone()));

                msets.setSpeciesFilepath(currentPath + "output" + File.separator + "maxent" + File.separator +
                        getName() + File.separator + "species_points.csv");
                msets.setOutputPath(currentPath + "output" + File.separator + "maxent" + File.separator + getName() +
                        File.separator);
                if (chkJackknife != null) {
                    msets.setDoJackknife(true);
                }
                if (chkResponseCurves != null) {
                    msets.setDoResponsecurves(true);
                }

                MaxentServiceImpl maxent = new MaxentServiceImpl();
                maxent.setMaxentSettings(msets);

                setStage(1);

                setProgress(0, "running Maxent");

                int exitValue = maxent.process(this);

                setProgress(1, "Maxent finished with exit value=" + exitValue);

                setStage(2);

                setProgress(0, "exporting results");

                Hashtable htProcess = new Hashtable();

                String[] imgExtensions = {".png", "_only.png", "_only_thumb.png", "_thumb.png"};

                if (isCancelled()) {
                    //
                } else if (exitValue == 0) {

                    // check if there is an error
                    String maxentError = getMaxentError(new File(msets.getOutputPath() + "maxent.log"), 2);
                    if (maxentError != null) {
                        setProgress(1, "failed: " + maxentError);
                        setCurrentState(FAILED);
                        if (maxentError.equals("Warning: Skipping species because it has 0 test samples")) {
                            setMessage("Warning: Skipping species because it has 0 test samples." +
                                    (msets.getRandomTestPercentage() > 0 ? "\nHint: Try to set the test percetage to '0'" : ""));
                        } else if (maxentError.equals("No species selected")) {
                            setMessage("No species selected.\nHint: Make sure your active area includes species occurrences");
                        }
                    } else {
                        // rename the env filenames to their display names
                        String pth_plots = currentPath + "output" + File.separator + "maxent" + File.separator +
                                getName() + File.separator + "plots" + File.separator;
                        String pth = currentPath + "output" + File.separator + "maxent" + File.separator + getName() +
                                File.separator;

                        readReplace(pth + "species.html", "Maxent model for species", "Maxent model for " + taxon);

                        String paramlist = "Model reference number: " + getName()
                                + "<br>Species: " + taxon
                                + "<br>Layers: <ul>";

                        LayerDAO layerDao = Client.getLayerDao();
                        for (int ei = 0; ei < envnameslist.length; ei++) {
                            System.out.println("LAYER NAME: " + envnameslist[ei]);
                            Layer lyr = layerDao.getLayerByName(envnameslist[ei]);

                            if (lyr != null) {
                                paramlist += "<li>" + lyr.getDisplayname() + " (" + envnameslist[ei] + ")</li>";
                            } else {
                                paramlist += "<li>" + displayNames[ei] + "</li>";
                            }

                            readReplace(pth + "species.html", "<td>" + envnameslist[ei] + "</td>", "<td>" + displayNames[ei] + "</td>");
                            readReplaceAfter(pth + "species.html", "(all continuous)", envnameslist[ei], displayNames[ei]);
                        }
                        paramlist += "</ul>";

                        readReplace(pth + "species.html", "end of this page.<br>", "end of this page.<br><p>" + paramlist + "</p>");
                        readReplace(pth + "species.html", msets.getOutputPath(), "");
                        readReplaceBetween(pth + "species.html", "Command line", "<br>", "");
                        readReplaceBetween(pth + "species.html", "Command line", "<br>", "");

                        // replace the summary
                        readReplace(pth + "species.html", "This page contains some analysis of the Maxent model for", "This <a href='http://www.cs.princeton.edu/~schapire/maxent/'>Maxent</a> v3.3.3e predictive model for");
                        readReplace(pth + "species.html", ", created", " was created");
                        readReplace(pth + "species.html", " using Maxent version 3.3.3e.", ".");
                        readReplace(pth + "species.html", "If you would like to do further analyses, the raw data used here is linked to at the end of this page", "Links at the bottom of this page to the raw data may be used for further analysis");

                        if (chkResponseCurves != null) {
                            StringBuffer sbTable = new StringBuffer();
                            String[] ctxlist = msets.getEnvVarToggler().split(" ");
                            if (msets.getEnvVarToggler().length() > 0) {
                                sbTable.append("<pre>");
                                for (String ctx : ctxlist) {
                                    if (ctx.startsWith("aloc_")) {
                                        continue;
                                    }
                                    sbTable.append("<span style='font-weight: bold; text-decoration: underline'>" + ctx + " legend</span><br />");
                                    sbTable.append(IOUtils.toString(new FileInputStream(GridCutter.getLayerPath(resolution, ctx) + ".txt")));
                                    sbTable.append("<br /><br />");
                                }
                                sbTable.append("</pre>");
                                readReplace(pth + "species.html", "<br><HR><H2>Analysis of variable contributions</H2><br>", sbTable.toString() + "<br><HR><H2>Analysis of variable contributions</H2><br>");
                            }
                        }

                        readReplaceBetween(pth + "species.html", "<br>Click <a href=species_explain.bat", "memory.<br>", "");
                        readReplaceBetween(pth + "species.html", "(A link to the Explain", "additive models.)", "");

                        StringBuffer removedSpecies = new StringBuffer();
                        try {
                            br = new BufferedReader(new FileReader(
                                    currentPath + "output" + File.separator + "maxent" + File.separator + getName() +
                                            File.separator + "Prediction_removedSpecies.txt"));
                            String ss;
                            while ((ss = br.readLine()) != null) {
                                removedSpecies.append(ss);
                            }
                            br.close();
                        } catch (Exception e) {
                        }
                        if (removedSpecies.length() > 0) {
                            String header = "'Sensitive species' have been masked out of the model. See: http://www.ala.org.au/about/program-of-projects/sds/\r\n\r\nLSID,Species scientific name,Taxon rank";
                            writeToFile(header + removedSpecies.toString(),
                                    currentPath + "output" + File.separator + "maxent" + File.separator + getName() + File.separator + "Prediction_maskedOutSensitiveSpecies.csv");

                            String insertBefore = "<a href = \"species.asc\">The";
                            String insertText = "<b><a href = \"Prediction_maskedOutSensitiveSpecies.csv\">'Sensitive species' masked out of the model</a></br></b>";
                            readReplace(pth + "species.html", insertBefore, insertText + insertBefore);
                        }

                        writeProjectionFile(msets.getOutputPath());

                        // if generated successfully, then add it to geoserver
                        String url = AlaspatialProperties.getGeoserverUrl() + "/rest/workspaces/ALA/coveragestores/maxent_" + getName() + "/file.arcgrid?coverageName=species_" + getName();
                        String extra = "";
                        String username = AlaspatialProperties.getGeoserverUsername();
                        String password = AlaspatialProperties.getGeoserverPassword();

                        // first zip up the file as it's going to be sent as binary
                        String[] infiles = {msets.getOutputPath() + "species.asc", msets.getOutputPath() + "species.prj"};
                        String ascZipFile = msets.getOutputPath() + "species.zip";
                        Zipper.zipFiles(infiles, ascZipFile);

                        // Upload the file to GeoServer using REST calls
                        UploadSpatialResource.loadResource(url, extra, username, password, ascZipFile);

                        htProcess.put("status", "success"); ///
                        htProcess.put("pid", getName());
                        htProcess.put("info", "/output/maxent/" + getName() + "/species.html");

                        //convert .asc to .grd/.gri
                        convertAscToDiva(msets.getOutputPath() + "species.asc", msets.getOutputPath() + getName());

                        setStage(3);

                        // generate the readme.txt file
                        CitationService.generatePredictionReadme(msets.getOutputPath(), msets.getSpeciesFilepath().substring(msets.getSpeciesFilepath().lastIndexOf("points")));

                        setProgress(1, "finished");

                        setCurrentState(SUCCESSFUL);

                        //write out infor for adjusting input parameters
                        System.out.println("MAXENT:" + cells + "," + envnameslist.length + " " + speciesCount + " " + (stageTimes[1] - stageTimes[0]) + " " + (stageTimes[2] - stageTimes[0]) + " " + (stageTimes[3] - stageTimes[2]));
                    }
                } else {
                    System.out.println("Failed 1");
                    setProgress(1, "failed");
                    setCurrentState(FAILED);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            setProgress(1, "failed: " + e.getMessage());
            setCurrentState(FAILED);
            setMessage("Error processing your Prediction request. Please try again or if problem persists, contact the Administrator.\n\nPlease quote the Prediction ID: " + getName());
        }
    }

    @Override
    public long getEstimate() {
        if (getProgress() == 0) {
            return 0;
        }

        long progTime;
        synchronized (progress) {
            progTime = progressTime;
        }
        long timeRemaining = 0;
        long t1 = 0, t2 = 0, t3 = 0;

        if (stage <= 0) { //data load; 0 to 0.2
            t1 += (cells * AlaspatialProperties.getAnalysisMaxentEstimateMult0()) * envnameslist.length; //default
            t1 = t1 + progTime - stageTimes[0];
        }
        if (stage <= 1) { //running; 0.2 to 0.9
            t2 += (cells * AlaspatialProperties.getAnalysisMaxentEstimateMult1()) * envnameslist.length; //default
            if (stage == 1) {
                t2 = t2 + progTime - stageTimes[1];
            }
        }
        if (stage > 1) { //data export + done
            t3 += 5000 * AlaspatialProperties.getAnalysisMaxentEstimateMult2(); //default
            if (stage == 2) {
                t3 = t3 + progTime - stageTimes[2];
            }
        }

        timeRemaining = t1 + t2 + t3;

        return timeRemaining;
    }

    @Override
    public double getProgress() {
        //return expected progress since cannot track internals

        long currentTime = System.currentTimeMillis();

        long progTime;
        synchronized (progress) {
            progTime = progressTime;
        }

        long t1 = 0, t2 = 0, t3 = 0;
        double d1, d2, d3;

        //progress is [time passed] / [time expected]
        if (stage <= 0) { //data load; 0 to 0.2
            t1 += (cells * AlaspatialProperties.getAnalysisMaxentEstimateMult0()) * envnameslist.length; //default
            d1 = (currentTime - stageTimes[0]) / (double) t1;
            if (d1 > 0.9) {
                d1 = 0.9;
            }
            d1 *= 0.2; //range limit
        } else {
            d1 = 0.2;
        }
        if (stage <= 1) { //running; 0.2 to 0.9
            t2 += (cells * AlaspatialProperties.getAnalysisMaxentEstimateMult1()) * envnameslist.length; //default
            if (stage == 1) {
                d2 = (currentTime - stageTimes[1]) / (double) t2;
            } else {
                d2 = 0;
            }
            if (d2 > 0.9) {
                d2 = 0.9;
            }
            d2 *= 0.7; //range limit
        } else {
            d2 = 0.7;
        }
        if (stage > 1) { //data export + done
            t3 += 5000 * AlaspatialProperties.getAnalysisMaxentEstimateMult2(); //default
            if (stage == 2) {
                d3 = (currentTime - stageTimes[2]) / (double) t3;
            } else {
                d3 = 0;
            }
            if (d3 > 0.9) {
                d3 = 0.9;
            }
            d3 *= 0.1; //range limit
        } else {
            d3 = 0.1;
        }

        return d1 + d2 + d3;
    }

    @Override
    public void setProgress(double d) {
        if (stage == 0) { //data load; 0 to 0.2
            progress = d / 5.0;
        } else if (stage == 1) { //running; 0.2 to 0.9
            progress = 0.2 + 10 * d / 7.0;
        } else { //exporting/done
            progress = 0.9 + d / 10.0;
        }
        super.setProgress(progress);
    }

    @Override
    public String getStatus() {
        if (getProgress() < 1) {
            String msg;
            if (stage == 0) { //data load; 0 to 0.2
                msg = "Data preparation, ";
            } else if (stage == 1) { //seeding; 0.2 to 0.9
                msg = "Running, ";
            } else {    //transforming data; 0.9 to 1.0
                msg = "Exporting results, ";
            }
            return msg + "est remaining: " + getEstimateInMinutes() + " min";
        } else {
            if (stage == -1) {
                return "not started, est: " + getEstimateInMinutes() + " min";
            } else {
                return "finished, total run time=" + Math.round(getRunTime() / 1000) + "s";
            }
        }
    }

    @Override
    public void setStage(int i) {
        super.setStage(i);
        if (i < 4) {
            stageTimes[i] = System.currentTimeMillis();
        }
    }

    public void setCells(int i) {
        cells = i;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(getName());
        sb.append("; Maxent");
        sb.append("; state=").append(getCurrentState());
        sb.append("; status=").append(getStatus());
        sb.append("; grid cell count=").append(cells);
        sb.append("; number of layers=").append(envnameslist.length);

        return sb.toString();
    }

    public void readReplaceBetween(String fname, String startOldText, String endOldText, String replText) {
        String line;
        StringBuffer sb = new StringBuffer();
        try {
            FileInputStream fis = new FileInputStream(fname);
            BufferedReader reader = new BufferedReader(new InputStreamReader(fis));
            while ((line = reader.readLine()) != null) {
                sb.append(line + "\n");
            }
            int start, end;
            start = sb.indexOf(startOldText);
            if (start >= 0) {
                end = sb.indexOf(endOldText, start + 1);
                sb.replace(start, end + endOldText.length(), replText);
            }
            reader.close();
            BufferedWriter out = new BufferedWriter(new FileWriter(fname));
            out.write(sb.toString());
            out.close();
        } catch (Throwable e) {
            System.err.println("*** exception ***");
            e.printStackTrace(System.out);
        }
    }

    private String setupSpecies(String speciesList, String outputpath) {
        try {
            File fDir = new File(outputpath);
            fDir.mkdir();

            File spFile = new File(fDir, "species_points.csv");
            PrintWriter spWriter = new PrintWriter(new BufferedWriter(new FileWriter(spFile)));

            spWriter.write(speciesList);
            spWriter.close();

            return spFile.getAbsolutePath();
        } catch (IOException ex) {
            ex.printStackTrace(System.out);
        }

        return null;
    }

    private void writeToFile(String text, String filename) {
        try {
            FileWriter fw = new FileWriter(filename);
            fw.append(text);
            fw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void writeProjectionFile(String outputpath) {
        try {
            File fDir = new File(outputpath);
            fDir.mkdir();

            PrintWriter spWriter = new PrintWriter(new BufferedWriter(new FileWriter(outputpath + "species.prj")));

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

            spWriter.write(sbProjection.toString());
            spWriter.close();

        } catch (IOException ex) {
            ex.printStackTrace(System.out);
        }
    }

    public String getImage() {
        return "output/maxent/" + getName() + "/plots/species_hidden.png";
    }

    AnalysisJob copy() {
        return new AnalysisJobMaxent(String.valueOf(System.currentTimeMillis()),
                currentPath, taxon, envlist, region, envelope, txtTestPercentage, chkJackknife, chkResponseCurves, resolution);
    }

    private String getMaxentError(File file, int count) {
        try {
            RandomAccessFile rf = new RandomAccessFile(file, "r");


            // first check if maxent threw a 'No species selected' error
            String nosp = rf.readLine(); // first line: date/time
            nosp = rf.readLine(); // second line: maxent version
            nosp = rf.readLine(); // third line: "No species selected"
            if (nosp.equals("No species selected")) {
                return "No species selected";
            }

            long flen = file.length() - 1;
            int nlcnt = -1;
            StringBuilder lines = new StringBuilder();
            while (nlcnt != count) {
                rf.seek(flen--);
                char c = (char) rf.read();
                lines.append(c);
                if (c == '\n') {
                    nlcnt++;
                }

            }
            String line = lines.reverse().toString();
            if (line.contains("Warning: Skipping species because it has 0 test samples")) {
                return "Warning: Skipping species because it has 0 test samples";
            }

            rf.close();
        } catch (Exception e) {
            System.out.println("Unable to read lines");
            e.printStackTrace(System.out);
        }

        // return false anyways
        return null;
    }

    private void convertAscToDiva(String asc, String grd) {
        try {
            //read asc
            BufferedReader br = new BufferedReader(new FileReader(asc));
            String s;

            //maxent output grid is:
            s = br.readLine();
            int ncols = Integer.parseInt(s.replace("ncols", "").trim());

            s = br.readLine();
            int nrows = Integer.parseInt(s.replace("nrows", "").trim());

            s = br.readLine();
            double lng1 = Double.parseDouble(s.replace("xllcorner", "").trim());

            s = br.readLine();
            double lat1 = Double.parseDouble(s.replace("yllcorner", "").trim());

            s = br.readLine();
            double div = Double.parseDouble(s.replace("cellsize", "").trim());

            s = br.readLine();
            double nodata = Double.parseDouble(s.replace("NODATA_value", "").trim());

            double[] data = new double[ncols * nrows];
            for (int i = 0; i < ncols * nrows; i++) {
                data[i] = Double.NaN;
            }
            int r = 0;
            while ((s = br.readLine()) != null) {
                String[] row = s.split(" ");
                for (int i = 0; i < row.length && i < ncols; i++) {
                    double v = Double.parseDouble(row[i]);
                    if (v != nodata) {
                        data[r * ncols + i] = v;
                    }
                }
                r++;
                if (r == nrows) {
                    break;
                }
            }
            br.close();

            Grid g = new Grid(null);
            g.writeGrid(grd, data, lng1, lat1, lng1 + ncols * div, lat1 + nrows * div, div, div, nrows, ncols);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private String writeFile(String contents, String outputpath, String filename) {
        try {
            File fDir = new File(outputpath);
            fDir.mkdir();

            File spFile = new File(fDir, filename);
            PrintWriter spWriter = new PrintWriter(new BufferedWriter(new FileWriter(spFile)));

            spWriter.write(contents);
            spWriter.close();

            return spFile.getAbsolutePath();
        } catch (IOException ex) {
            ex.printStackTrace(System.out);
        }

        return null;
    }
}
