/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ala.spatial.web.zk;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import org.ala.spatial.util.TabulationSettings;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.springframework.util.StringUtils;
import org.zkoss.zk.ui.Component;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.util.GenericForwardComposer;
import org.zkoss.zul.Button;
import org.zkoss.zul.Grid;
import org.zkoss.zul.Label;
import org.zkoss.zul.ListModelList;
import org.zkoss.zul.Row;
import org.zkoss.zul.RowRenderer;

/**
 *
 * @author ajay
 */
public class LogsZK extends GenericForwardComposer {

    private Label lblMessage;
    private Button refresh;
    private Grid grid;
    final private static String LOGS_BASE_DIR = "/data/logs/";
    ArrayList logList = new ArrayList();
    private static List processedfiles = null;
    String[] headers = {"date", "servername", "userip", "useremail", "processid", "sessionid", "actiontype", "lsid", "layers", "method", "params", "downloadfile", "message"};
    static String ALA_LOGGING_SERVICE_URL = "http://ala-biocache1.vm.csiro.au:8080/ala-logger-service/service/logger/";

    @Override
    public void doAfterCompose(Component comp) throws Exception {
        super.doAfterCompose(comp);

        //lblMessage.setValue("Log Analyser");
        System.out.println("Calling loadLogs");

        //grid.setVflex(true);
        //grid.setMold("paging");        

        loadLogs();
    }

    public void onClick$refresh(Event e) {
        System.out.println("Refreshing logs");
        //grid.setMold(null);
        //grid.getChildren().clear();
        logList.clear();
        //grid.setMold("paging");
        loadLogs();

    }

    private void loadLogs() {
        try {
            File logsDir = new File(LOGS_BASE_DIR); // /Library/Tomcat/Home/logs/

            FileFilter ff = new WildcardFileFilter("useractions.*"); // useractions.log.2011-01-1*
            File[] files = logsDir.listFiles(ff);
            System.out.println("Data...");
            for (int i = 0; i < files.length; i++) {
                System.out.println(files[i]);

                CSVReader reader = new CSVReader(new FileReader(files[i]));
                logList.addAll(reader.readAll());

            }

            ListModelList lm = new ListModelList(logList);
            grid.setModel(lm);
            grid.setRowRenderer(new RowRenderer() {

                @Override
                public void render(Row row, Object data) throws Exception {
                    String[] recs = (String[]) data;
                    for (int i = 0; i < recs.length; i++) {
                        String col = recs[i].trim();
                        if (col.startsWith("\"")) {
                            col = col.substring(1);
                        }
                        new Label(col).setParent(row);
                    }

                }
            });


        } catch (Exception e) {
            System.out.println("Error loading logs: ");
            e.printStackTrace(System.out);
        }
    }

    private static void postInfoToLogger(String ip, Hashtable<String, Integer> uidCounts) {
        try {

            StringBuilder sbInfo = new StringBuilder();
            sbInfo.append("{");
            sbInfo.append("\"eventTypeId\": 1002,");
            sbInfo.append("\"comment\": \"Sampling download from ALA Spatial Portal\",");
            sbInfo.append("\"userEmail\" : \"spatial@ala\",");
            sbInfo.append("\"userIP\": \"" + ip + "\",");
            sbInfo.append("\"recordCounts\" : {");
            Iterator<String> uit = uidCounts.keySet().iterator();
            while (uit.hasNext()) {
                String uid = uit.next();
                int cnt = uidCounts.get(uid).intValue();
                //System.out.println(uid + " - " + cnt);

                sbInfo.append("\"" + uid + "\":" + cnt);
                if (uit.hasNext()) {
                    sbInfo.append(",");
                }
            }

            sbInfo.append("}");
            sbInfo.append("}");

            System.out.println("Sending: \n" + sbInfo.toString());

            HttpClient client = new HttpClient();

            // http://ala-biocache1.vm.csiro.au:8080/ala-logger-service/service/logger/
            PostMethod post = new PostMethod(ALA_LOGGING_SERVICE_URL);

            RequestEntity entity = new StringRequestEntity(sbInfo.toString(), "application/json", "utf-8");
            post.setRequestEntity(entity);


            System.out.println("Sending to ALA-Logger...");
            int result = client.executeMethod(post);

            writeLoggedInfo(sbInfo.toString());

            if (result == HttpStatus.SC_OK) {
                System.out.println("Successfully sent logging service:");
                System.out.println(post.getResponseBodyAsString());
            } else {
                System.out.println("Not a clean response: " + result);
            }
        } catch (Exception ex) {
            System.out.println("postInfo.error:");
            ex.printStackTrace(System.out);
        }
    }

    private static void readProcessedFile() {
        try {
            File procFile = new File(LOGS_BASE_DIR + "useraction.processed");
            if (!procFile.exists()) {
                procFile.createNewFile();
            }
            CSVReader reader = new CSVReader(new FileReader(procFile));
            processedfiles = reader.readAll();
            reader.close();
        } catch (Exception e) {
            System.out.println("Error reading processed file");
            e.printStackTrace(System.out);
        }
    }

    private static void writeProcessedFile() {
        try {
            File procFile = new File(LOGS_BASE_DIR + "useraction.processed");
            CSVWriter writer = new CSVWriter(new FileWriter(procFile));
            writer.writeAll(processedfiles);
            writer.close();
        } catch (Exception e) {
        }
    }

    private static void writeLoggedInfo(String info) {
        try {
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(LOGS_BASE_DIR + "useraction.sent", true)));
            out.write(info);
            out.close();
        } catch (Exception e) {
            System.out.println("Error, unable to write useraction.sent");
            e.printStackTrace(System.out);
        }
    }

    public static void processLogFiles() {
        processLogFiles(ALA_LOGGING_SERVICE_URL);
    }

    public static void processLogFiles(String url) {
        try {

            ALA_LOGGING_SERVICE_URL = url;

            readProcessedFile();

            File logsDir = new File(LOGS_BASE_DIR); // /Library/Tomcat/Home/logs/

            FileFilter ff = new WildcardFileFilter("useractions.dev.log"); // useractions.log.2011-01-1*
            File[] files = logsDir.listFiles(ff);
            System.out.println("Data...");
            for (int i = 0; i < files.length; i++) {
                System.out.println(files[i]);

                //if (processedfiles.contains(files[i].getAbsolutePath())) continue;
                Iterator pit = processedfiles.iterator();
                boolean isprocessed = false;
                while (pit.hasNext()) {
                    String[] rows = (String[]) pit.next();
                    if (rows[0].contains(files[i].getAbsolutePath())) {
                        isprocessed = true;
                        break;
                    }
                }

                if (isprocessed) {
                    continue;
                }

                CSVReader reader = new CSVReader(new FileReader(files[i]));
                List rows = reader.readAll();
                Iterator it = rows.iterator();

                while (it.hasNext()) {
                    String[] row = (String[]) it.next();
                    String ip = (row[2].trim().startsWith("\"")) ? row[2].trim().substring(1) : row[2].trim();
                    String method = (row[9].trim().startsWith("\"")) ? row[9].trim().substring(1) : row[9].trim();
                    String filename = row[11].trim();
                    if (filename.endsWith(".zip")) {
                        filename = filename.trim().substring(filename.lastIndexOf("/") + 1);
                    }

                    System.out.println("Method: " + method);
                    System.out.println("filename: " + filename);

                    if (method.equals("Sampling") && filename.endsWith(".zip")) {
                        System.out.println("Got a record, looking for " + filename);
                        String zipfilepath = TabulationSettings.base_output_dir + "output/sampling/" + filename;
                        System.out.println("Details for " + ip + " from " + zipfilepath);
                        File zipfile = new File(zipfilepath);
                        if (zipfile.exists() && (zipfile.length() > 0)) {
                            ZipFile zf = new ZipFile(zipfile);
                            Hashtable<String, Integer> uidCounts = readZipFile(zf);
                            postInfoToLogger(ip, uidCounts);
                        }
                    } else if (method.equals("Sampling") && filename.trim().equals("")) {
                        System.out.println("Only sampling, no file available");
                    } else {
                        System.out.println("no sampling and zip file");
                    }

                }

                reader.close();

                String[] fileproc = {files[i].getAbsolutePath()};
                processedfiles.add(fileproc);
            }

            writeProcessedFile();

        } catch (Exception e) {
            System.out.println("Error loading logs: ");
            e.printStackTrace(System.out);
        }
    }

    public static void processLogFilesAsLines(String url) {
        try {
            ALA_LOGGING_SERVICE_URL = url;

            readProcessedFile();

            File logsDir = new File(LOGS_BASE_DIR); // /Library/Tomcat/Home/logs/

            FileFilter ff = new WildcardFileFilter("useractions.dev.log"); // useractions.log.2011-01-1*
            File[] files = logsDir.listFiles(ff);
            System.out.println("Data...");
            for (int i = 0; i < files.length; i++) {

                System.out.println(files[i]);

                //if (processedfiles.contains(files[i].getAbsolutePath())) continue;
                Iterator pit = processedfiles.iterator();
                boolean isprocessed = false;
                while (pit.hasNext()) {
                    String[] rows = (String[]) pit.next();
                    if (rows[0].contains(files[i].getAbsolutePath())) {
                        isprocessed = true;
                        break;
                    }
                }

                if (isprocessed) {
                    continue;
                }

                //CSVReader reader = new CSVReader(new FileReader(files[i]));
                BufferedReader reader = new BufferedReader(new FileReader(files[i]));
                String line = null;
                while ((line = reader.readLine()) != null) {
                    //System.out.println(line);

                    if (line.contains("Sampling") && line.contains(".zip") && line.contains("output/sampling")) {
                        //System.out.println("\nFound sampling");
                        String[] fields = StringUtils.commaDelimitedListToStringArray(line);
                        String ip = (fields[3].trim().startsWith("\"")) ? fields[3].trim().substring(1) : fields[3].trim();
                        ip = ip.trim().replaceAll("\"", ""); 
                        String filename = fields[fields.length - 2];
                        if (filename.contains(".zip")) {
                            filename = filename.substring(filename.lastIndexOf("/")+1).trim().replaceAll("\"", "");
                        }

                        System.out.println("Looking up " + filename + " for " + ip);
                        String zipfilepath = TabulationSettings.base_output_dir + "output/sampling/" + filename;
                        System.out.println("Details for " + ip + " from " + zipfilepath);
                        File zipfile = new File(zipfilepath);
                        if (zipfile.exists() && (zipfile.length() > 0)) {
                            ZipFile zf = new ZipFile(zipfile);
                            Hashtable<String, Integer> uidCounts = readZipFile(zf);
                            postInfoToLogger(ip, uidCounts);
                        }

                    }

                    //System.out.println("\n\n");
                }


                String[] fileproc = {files[i].getAbsolutePath()};
                processedfiles.add(fileproc);
            }

            writeProcessedFile();

        } catch (Exception e) {
        }

    }

    private static void processSamplesDownloads() {
        try {
            File downloadDir = new File(TabulationSettings.base_output_dir + "sampling/");
            FileFilter ff = new WildcardFileFilter("*.zip");
            File[] files = downloadDir.listFiles(ff);

            for (int i = 0; i < files.length; i++) {
                System.out.println("Loading " + files[i]);
                ZipFile zf = new ZipFile(files[i]);
                readZipFile(zf);
            }

        } catch (Exception e) {
            System.out.println("Exception processing samples folder");
            e.printStackTrace(System.out);
        }
    }

    private static Hashtable<String, Integer> readZipFile(ZipFile zf) {
        try {

            Hashtable<String, Integer> uidCounts = new Hashtable<String, Integer>();

            Enumeration entries = zf.entries();
            while (entries.hasMoreElements()) {
                ZipEntry ze = (ZipEntry) entries.nextElement();

                if (ze.getName().contains("samples.csv")) {
                    BufferedReader br = new BufferedReader(
                            new InputStreamReader(zf.getInputStream(ze)));

                    CSVReader reader = new CSVReader(br);
                    List rows = reader.readAll();
                    Iterator it = rows.iterator();

                    // get the header and look for various _uid columns
                    Vector uidColumns = new Vector();
                    if (it.hasNext()) {
                        String[] header = (String[]) it.next();
                        for (int i = 0; i < header.length; i++) {
                            if (header[i].endsWith("_uid")) {
                                uidColumns.add(i);
                            }
                        }
                    }

                    for (; it.hasNext();) {
                        String[] row = (String[]) it.next();

                        Enumeration uids = uidColumns.elements();
                        while (uids.hasMoreElements()) {
                            int uid = ((Integer) uids.nextElement()).intValue();

                            int uc = 0;
                            if (uidCounts.containsKey(row[uid])) {
                                uc = uidCounts.get(row[uid]).intValue();
                                //System.out.println(row[uid] + " available: " + uc);
                            } else {
                                //System.out.println(row[uid] + " is new");
                            }
                            uc++;
                            if (!row[uid].trim().equals("")) {
                                uidCounts.put(row[uid], new Integer(uc));
                            }

                            uc = 0;
                        }
                    }

                    break;
                }
            }

            return uidCounts;

        } catch (Exception e) {
            System.out.println("Exception in readZipFile: ");
            e.printStackTrace(System.out);
        }

        return null;
    }

    public static void main(String[] args) {
        processLogFilesAsLines(ALA_LOGGING_SERVICE_URL);
    }
}
