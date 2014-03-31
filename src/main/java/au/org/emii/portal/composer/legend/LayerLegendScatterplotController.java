/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package au.org.emii.portal.composer.legend;

import au.org.ala.spatial.data.Query;
import au.org.ala.spatial.data.ScatterplotData;
import au.org.ala.spatial.sampling.Sampling;
import au.org.ala.spatial.util.CommonData;
import au.org.ala.spatial.util.LayersUtil;
import au.org.ala.spatial.util.SelectedArea;
import au.org.ala.spatial.util.Util;
import au.org.emii.portal.composer.MapComposer;
import au.org.emii.portal.composer.UtilityComposer;
import au.org.emii.portal.menu.MapLayer;
import au.org.emii.portal.settings.SettingsSupplementary;
import au.org.emii.portal.util.LayerUtilities;
import org.ala.layers.legend.Facet;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.codehaus.jackson.map.ObjectMapper;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.util.Clients;
import org.zkoss.zul.*;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Adam
 */
public class LayerLegendScatterplotController extends UtilityComposer {

    private static Logger logger = Logger.getLogger(LayerLegendScatterplotController.class);

    private static final String NUMBER_SERIES = "Number series";
    private static final String ACTIVE_AREA_SERIES = "In Active Area";
    private SettingsSupplementary settingsSupplementary = null;

    Textbox tbxChartSelection;
    Label tbxSelectionCount;
    Label tbxRange;
    Label tbxDomain;
    Label tbxMissingCount;
    ScatterplotData data;
    LayersUtil layersUtil;
    Div scatterplotButtons;
    Div scatterplotDownloads;
    Div divHighlightArea;
    MapLayer mapLayer = null;
    Checkbox chkSelectMissingRecords;

    Label lblMissing;
    Button addNewLayers;
    Combobox cbHighlightArea;

    ScatterplotLayerLegendComposer layerWindow = null;

    @Override
    public void afterCompose() {
        super.afterCompose();

        this.addEventListener("onSize", new EventListener() {

            @Override
            public void onEvent(Event event) throws Exception {
                redraw();
            }
        });
    }

    @Override
    public void doEmbedded() {
        super.doEmbedded();
        redraw();
    }

    @Override
    public void doOverlapped() {
        super.doOverlapped();
        redraw();
    }

    public ScatterplotData getScatterplotData() {
        if (data == null) {
            if (mapLayer == null) {
                data = new ScatterplotData();
            } else {
                data = mapLayer.getScatterplotData();
            }
        }
        return data;
    }

    public void onChange$tbxChartSelection(Event event) {
        try {
            //input order is [x, y, width, height]
            // + optional bounding box coordinates [x1,y1,x2,y2]
            logger.debug(event.getData());
            String[] coordsStr = ((String) event.getData()).replace("px", "").split(",");
            double[] coordsDbl = new double[coordsStr.length];
            for (int i = 0; i < coordsStr.length; i++) {
                coordsDbl[i] = Double.parseDouble(coordsStr[i]);
            }

            String params = "?minx=" + coordsDbl[0] + "&miny=" + coordsDbl[1] + "&maxx=" + coordsDbl[2] + "&maxy=" + coordsDbl[3];

            data.imagePath = "http://localhost:8082/alaspatial/ws/scatterplot/img/" + data.getId() + params;

            ObjectMapper om = new ObjectMapper();
            Map map = om.readValue(new URL(data.imagePath), Map.class);

            data.prevSelection = new double[4];
            data.prevSelection[0] = Double.parseDouble((String) map.get("minx"));
            data.prevSelection[1] = Double.parseDouble((String) map.get("miny"));
            data.prevSelection[2] = Double.parseDouble((String) map.get("maxx"));
            data.prevSelection[3] = Double.parseDouble((String) map.get("maxy"));

            Facet f = getFacetIn();
            if (f != null) {
                mapLayer.setHighlight(f.toString());
            } else {
                mapLayer.setHighlight(null);
            }

            getMapComposer().applyChange(mapLayer);

            tbxChartSelection.setText("");
            tbxDomain.setValue(String.format("%s: %g - %g", data.getLayer1Name(), data.prevSelection[1], data.prevSelection[3]));
            tbxRange.setValue(String.format("%s: %g - %g", data.getLayer2Name(), data.prevSelection[0], data.prevSelection[2]));

            data.imagePath = null;
            redraw();
        } catch (Exception e) {
            logger.error("failed to build scatterplot legend", e);
            clearSelection();
            getMapComposer().applyChange(mapLayer);
        }
    }

    void redraw() {
        int width = Integer.parseInt(this.getWidth().replace("px", "")) - 20;
        int height = Integer.parseInt(this.getHeight().replace("px", "")) - Integer.parseInt(tbxChartSelection.getHeight().replace("px", ""));
        if (height > width) {
            height = width;
        } else {
            width = height;
        }

        String script = "updateScatterplot(" + width + "," + height + ",'url(" + data.imagePath + ")')";
        Clients.evalJavaScript(script);

        scatterplotDownloads.setVisible(true);

        if (data.getMissingCount() > 0) {
            tbxMissingCount.setValue("(" + data.getMissingCount() + ")");
            chkSelectMissingRecords.setVisible(true);
        } else {
            tbxMissingCount.setValue("");
            chkSelectMissingRecords.setVisible(false);
        }
    }

    private void registerScatterPlotSelection() {
        try {
            double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
            if (data.prevSelection != null) {
                x1 = data.prevSelection[0];
                x2 = data.prevSelection[1];
                y1 = data.prevSelection[2];
                y2 = data.prevSelection[3];
            }

            if (data.getLayer1() != null && data.getLayer1().length() > 0
                    && data.getLayer2() != null && data.getLayer2().length() > 0) {
                Facet f = getFacetIn();

                int count = 0;
                if (f != null) {
                    Query q = data.getQuery().newFacet(f, false);
                    count = q.getOccurrenceCount();
                }
                updateCount(String.valueOf(count));
            }
        } catch (Exception e) {
            logger.error("error updating scatterplot selection", e);
            clearSelection();
            getMapComposer().applyChange(mapLayer);
        }
    }

    void updateCount(String txt) {
        try {
            data.selectionCount = Integer.parseInt(txt);
            tbxSelectionCount.setValue("Records selected: " + txt);
            if (data.selectionCount > 0) {
                addNewLayers.setVisible(true);
            } else {
                addNewLayers.setVisible(false);
            }
            scatterplotButtons.setVisible(true);
        } catch (Exception e) {
        }
    }

    void clearSelection() {
        tbxSelectionCount.setValue("");
        addNewLayers.setVisible(false);
        tbxRange.setValue("");
        tbxDomain.setValue("");

        data.prevSelection = null;

        chkSelectMissingRecords.setChecked(false);

        getScatterplotData().setEnabled(false);

        scatterplotDownloads.setVisible(false);

        data.prevSelection = null;

        scatterplotButtons.setVisible(false);
    }

    public void onClick$addSelectedRecords(Event event) {
        Facet f = getFacetIn();
        if (f != null) {
            addUserLayer(data.getQuery().newFacet(getFacetIn(), true), "IN " + data.getSpeciesName(), "from scatterplot in group", data.selectionCount);
        }
    }

    public void onClick$addUnSelectedRecords(Event event) {
        addUserLayer(data.getQuery().newFacet(getFacetOut(), true), "OUT " + data.getSpeciesName(), "from scatterplot out group", data.results.split("\n").length - data.selectionCount - 1);   //-1 for header
    }

    void addUserLayer(Query query, String layername, String description, int numRecords) {
        layername = StringUtils.capitalize(layername);

        getMapComposer().mapSpecies(query, layername, "species", -1, LayerUtilities.SPECIES, null, -1, MapComposer.DEFAULT_POINT_SIZE,
                MapComposer.DEFAULT_POINT_OPACITY, Util.nextColour());
    }

    public void onClick$addNewLayers(Event event) {
        onClick$addUnSelectedRecords(null);
        onClick$addSelectedRecords(null);
    }

    public void onClick$scatterplotImageDownload(Event event) {
        try {
            Filedownload.save(new File(data.imagePath), "image/png");
        } catch (Exception e) {
            logger.error("error saving scatterplot image: " + data.imagePath, e);
        }
    }

    public void onClick$scatterplotDataDownload(Event event) {

        try {
            ObjectMapper om = new ObjectMapper();

            String csv = om.readValue(new URL("http://localhost:8082/alaspatial/ws/scatterplot/csv/" + data.getId()), String.class);

            Filedownload.save(csv, "text/plain", "scatterplot.csv");
        } catch (Exception e) {
            logger.error("error downloading scatterplot csv for id: " + data.getId(), e);
        }
    }

    public void onCheck$chkSelectMissingRecords(Event event) {
        try {
            registerScatterPlotSelection();

            ScatterplotData d = getScatterplotData();
            d.setEnabled(true);

            Facet f = getFacetIn();
            if (f == null) {
                mapLayer.setHighlight(null);
            } else {
                mapLayer.setHighlight(f.toString());
            }

            getMapComposer().applyChange(mapLayer);

            tbxChartSelection.setText("");

            data.imagePath = null;
            redraw();
        } catch (Exception e) {
            e.printStackTrace();
            clearSelection();
            getMapComposer().applyChange(mapLayer);
        }
    }

    public void updateFromLegend() {
        updateFromLegend(
                layerWindow.getRed(),
                layerWindow.getGreen(),
                layerWindow.getBlue(),
                layerWindow.getOpacity(),
                layerWindow.getPlotSize(),
                layerWindow.getColourMode());

        if (mapLayer != null) {
            mapLayer.setColourMode(layerWindow.getColourMode());
            mapLayer.setRedVal(layerWindow.getRed());
            mapLayer.setGreenVal(layerWindow.getGreen());
            mapLayer.setBlueVal(layerWindow.getBlue());
            mapLayer.setOpacity(layerWindow.getOpacity() / 100.0f);
            mapLayer.setSizeVal(layerWindow.getSize());

            getMapComposer().applyChange(mapLayer);
        }
    }

    public void updateFromLegend(int red, int green, int blue, int opacity, int size, String colourMode) {
        data.red = red;
        data.green = green;
        data.blue = blue;
        data.opacity = opacity;
        data.size = size;
        data.colourMode = colourMode;

        data.imagePath = null;
        redraw();
    }

    private void store() {
        if (mapLayer != null) {
            mapLayer.setScatterplotData(data);
        }
    }

    private Facet getFacetIn() {
        String fq = null;
        String e1 = CommonData.getLayerFacetName(data.getLayer1());
        String e2 = CommonData.getLayerFacetName(data.getLayer2());

        if (chkSelectMissingRecords.isChecked() && data.prevSelection == null) {
            fq = "-(" + e1 + ":[* TO *] AND " + e2 + ":[* TO *])";
        } else if (data.prevSelection != null) {
            double x1 = data.prevSelection[0];
            double x2 = data.prevSelection[1];
            double y1 = data.prevSelection[2];
            double y2 = data.prevSelection[3];

            Facet f1 = new Facet(e1, y1, y2, true);
            Facet f2 = new Facet(e2, x1, x2, true);

            if (chkSelectMissingRecords.isChecked()) {
                fq = "-(-(" + f1.toString() + " AND " + f2.toString() + ") AND " + e1 + ":[* TO *] AND " + e2 + ":[* TO *])";
            } else {
                fq = f1.toString() + " AND " + f2.toString();
            }
        }

        return Facet.parseFacet(fq);
    }

    private Facet getFacetOut() {
        String fq = "*:*";
        String e1 = CommonData.getLayerFacetName(data.getLayer1());
        String e2 = CommonData.getLayerFacetName(data.getLayer2());
        if (chkSelectMissingRecords.isChecked() && data.prevSelection == null) {
            fq = e1 + ":[* TO *] AND " + e2 + ":[* TO *]";
        } else if (data.prevSelection != null) {
            double x1 = data.prevSelection[0];
            double x2 = data.prevSelection[1];
            double y1 = data.prevSelection[2];
            double y2 = data.prevSelection[3];

            Facet f1 = new Facet(e1, y1, y2, true);
            Facet f2 = new Facet(e2, x1, x2, true);
            if (chkSelectMissingRecords.isChecked()) {
                fq = "-(" + f1.toString() + " AND " + f2.toString() + ") AND " + e1 + ":[* TO *] AND " + e2 + ":[* TO *]";
            } else {
                fq = "-(" + f1.toString() + " AND " + f2.toString() + ")";
            }
        }

        return Facet.parseFacet(fq);
    }

    private double[][] sample(double[] points) {
        double[][] p = new double[points.length / 2][2];
        for (int i = 0; i < points.length; i += 2) {
            p[i / 2][0] = points[i];
            p[i / 2][1] = points[i + 1];
        }

        ArrayList<String> layers = new ArrayList<String>();
        layers.add(CommonData.getLayerFacetName(data.getLayer1()));
        layers.add(CommonData.getLayerFacetName(data.getLayer2()));
        List<String[]> sample = Sampling.sampling(layers, p);

        double[][] d = new double[p.length][2];
        for (int j = 0; j < p.length; j++) {
            for (int i = 0; i < sample.size(); i++) {
                try {
                    d[j][i] = Double.parseDouble(sample.get(i)[j]);
                } catch (Exception e) {
                    d[j][i] = Double.NaN;
                }
            }
        }

        return d;
    }

    private void updateCbHighlightArea() {
        for (int i = cbHighlightArea.getItemCount() - 1; i >= 0; i--) {
            cbHighlightArea.removeItemAt(i);
        }

        boolean selectionSuccessful = false;
        for (MapLayer ml : getMapComposer().getPolygonLayers()) {
            Comboitem ci = new Comboitem(ml.getDisplayName());
            ci.setValue(ml);
            ci.setParent(cbHighlightArea);
            if (data != null && data.getHighlightSa() != null
                    && data.getHighlightSa().getMapLayer().getName().equals(ml.getName())) {
                cbHighlightArea.setSelectedItem(ci);
                selectionSuccessful = true;
            }
        }

        //this may be a deleted layer or current view or au or world
        if (!selectionSuccessful && data != null
                && data.getHighlightSa() != null) {
            MapLayer ml = data.getHighlightSa().getMapLayer();
            if (ml != null) {
                Comboitem ci = new Comboitem(ml.getDisplayName() + " (DELETED LAYER)");
                ci.setValue(ml);
                ci.setParent(cbHighlightArea);
                cbHighlightArea.setSelectedItem(ci);
            } else {
                String name = "Previous area";
                if (data.getHighlightSa().getWkt() != null) {
                    if (data.getHighlightSa().getWkt().equals(CommonData.AUSTRALIA_WKT)) {
                        name = "Australia";
                    } else if (data.getHighlightSa().getWkt().equals(CommonData.WORLD_WKT)) {
                        name = "World";
                    }
                }
                Comboitem ci = new Comboitem(name);
                ci.setValue(data.getHighlightSa().getWkt());
                ci.setParent(cbHighlightArea);
                cbHighlightArea.setSelectedItem(ci);
            }
        }
    }

    public void onSelect$cbHighlightArea(Event event) {
        if (cbHighlightArea.getSelectedItem() != null) {
            if (cbHighlightArea.getSelectedItem().getValue() instanceof MapLayer) {
                MapLayer ml = ((MapLayer) cbHighlightArea.getSelectedItem().getValue());
                SelectedArea sa = new SelectedArea(ml, ml.getWKT());
                data.setHighlightSa(sa);
            } else {
                String wkt = (String) cbHighlightArea.getSelectedItem().getValue();
                SelectedArea sa = new SelectedArea(null, wkt);
                data.setHighlightSa(sa);
            }
        } else {
            data.setHighlightSa(null);
        }
        data.imagePath = null;
        redraw();
    }

    public void onClick$bClearHighlightArea(Event event) {
        cbHighlightArea.setSelectedIndex(-1);
        onSelect$cbHighlightArea(null);
    }
}