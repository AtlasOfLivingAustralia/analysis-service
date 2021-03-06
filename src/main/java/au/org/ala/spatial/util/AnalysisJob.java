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

import org.springframework.web.client.RestOperations;
import org.springframework.web.client.RestTemplate;

import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import java.util.Map;

/**
 * Parent class for analysis jobs.
 *
 * @author Adam
 */
public class AnalysisJob extends Thread implements Serializable {
    public final static String WAITING = "WAITING";
    public final static String RUNNING = "RUNNING";
    public final static String SUCCESSFUL = "SUCCESSFUL";
    public final static String FAILED = "FAILED";
    public final static String CANCELLED = "CANCELLED";
    private final String DATE_FORMAT_NOW = "dd-MM-yyyy HH:mm:ss";
    Double progress;            //progress 0 to 1
    long progressTime;          //time of last set progress
    int stage;                  //for use by extended classes
    String status;              //status of thread
    StringBuffer log;           //log of job events
    long runTime;               //total run time in ms
    String currentState;        //state of job; WAITING, RUNNING, SUCCESSFUL, FAILED
    long[] estimatePairs = null;//pairs of time requested of time remaining
    String inputs;              //input name:value; pairs
    String message;

    public AnalysisJob(String pid) {
        this.setName(pid);
        status = "";
        log = new StringBuffer();
        stage = -1;
        progress = new Double(0);
        message = "";
        setCurrentState(WAITING);
    }

    public static void main(String[] args) {
        List<Map<String, String>> results = null;
        final String url = "http://bie.ala.org.au/ws/species/guids/bulklookup.json";
        try {
            //String jsonString="";
            String[] list = new String[]{
                    "urn:lsid:biodiversity.org.au:apni.taxon:306711",
                    "Grallina cyanoleuca|urn:lsid:biodiversity.org.au:afd.taxon:8f51b1e4-44c3-4a4c-9a25-a54d6a345ceb|Magpie-lark|ANIMALIA|MONARCHIDAE",
                    "||||"};
            RestOperations restTemplate = new RestTemplate();
            Map searchDTOList = restTemplate.postForObject(url, list, Map.class);
            //System.out.println(test);
            results = (List<Map<String, String>>) searchDTOList.get("searchDTOList");
        } catch (Exception ex) {
        }

    }

    public long smoothEstimate(long nextEstimate) {
        int estimates_length = (int) (AlaspatialProperties.getAnalysisEstimateSmoothing() * 2);
        long difference = -1 * progressTime;
        long currentTime = System.currentTimeMillis();
        difference += currentTime;

        if (estimatePairs == null) {
            estimatePairs = new long[estimates_length];
            for (int i = 0; i < estimates_length; i += 2) {
                estimatePairs[i] = currentTime;
                estimatePairs[i + 1] = nextEstimate;
            }
        }

        //shift estimates back one
        for (int i = 2; i < estimates_length; i += 2) {
            estimatePairs[i] = estimatePairs[i - 2];
            estimatePairs[i + 1] = estimatePairs[i - 1];
        }

        //apply this estimate
        estimatePairs[0] = currentTime;
        estimatePairs[1] = nextEstimate;

        //calculate
        double smoothEstimate = 0;
        for (int i = 0; i < estimates_length; i += 2) {
            long timeSinceEstimate = currentTime - estimatePairs[i];
            smoothEstimate += estimatePairs[i + 1] - timeSinceEstimate;
        }

        //return average time remaining
        long estimate = ((long) Math.round(smoothEstimate / (estimates_length / 2))) - difference;
        if (estimate < 1000) estimate = 1000;
        return estimate;
    }

    public long getEstimate() {
        return smoothEstimate(0);
    }

    public boolean isFinished() {
        return getCurrentState().equals(CANCELLED) ||
                getCurrentState().equals(SUCCESSFUL) ||
                getCurrentState().equals(FAILED);
    }

    public boolean isCancelled() {
        return getCurrentState().equals(CANCELLED);
    }

    public String getEstimateInMinutes() {
        double minutes = getEstimate() / 60000.0;
        if (minutes < 0.5) {
            return "< 1";
        } else {
            return String.valueOf((int) Math.ceil(minutes));
        }
    }

    public boolean cancel() {
        setStatus(CANCELLED);
        return true;
    }

    public String getLog() {
        return log.toString();
    }

    public void log(String s) {
        log.append("\r\n").append(now()).append("> ").append(s);
    }

    public long getRunTime() {
        return runTime;
    }

    public void setRunTime(long l) {
        runTime = l;
    }

    public String getStatus() {
        return status;
    }

    public void setStatus(String s) {
        status = s;
    }

    public void setProgress(double d, String s) {
        setProgress(d);
        log(s);
    }

    public double getProgress() {
        return progress;
    }

    public void setProgress(double d) {
        synchronized (progress) {
            progress = d;
            progressTime = System.currentTimeMillis();
        }
    }

    public int getStage() {
        return stage;
    }

    public void setStage(int i) {
        stage = i;
        setProgress(0); //for this stage
    }

    private String now() {
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
        return sdf.format(cal.getTime());
    }

    public String getCurrentState() {
        return currentState;
    }

    public void setCurrentState(String state) {
        currentState = state;
    }

    public String getInputs() {
        return inputs;
    }

    public void setInputs(String inputs_) {
        inputs = inputs_;
    }

    public String getMessage() {
        return message;
    }

    public void setMessage(String message) {
        this.message = message;
    }

    public String getImage() {
        return null;
    }

    AnalysisJob copy() {
        return new AnalysisJob(String.valueOf(System.currentTimeMillis()));
    }
}
