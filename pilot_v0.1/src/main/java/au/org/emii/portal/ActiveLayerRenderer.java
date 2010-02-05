package au.org.emii.portal;

import org.zkoss.zk.ui.Executions;
import org.zkoss.zul.Checkbox;
import org.zkoss.zul.Image;
import org.zkoss.zul.Label;
import org.zkoss.zul.Listcell;
import org.zkoss.zul.Listitem;
import org.zkoss.zul.ListitemRenderer;
import org.zkoss.zul.Popup;


public class ActiveLayerRenderer implements ListitemRenderer {

		
	public void render(Listitem item, Object data) throws Exception {
		final MapLayer layer = (MapLayer) data;
		Listcell listcell = new Listcell();
		Checkbox checkbox = new Checkbox();
		
		/*
		 * In the past it was assumed that we just set true here - this is not the
		 * case because this method is called to re-render the list after the 
		 * user has been fiddling with the checkboxes in some cases (dnd events) 
		 */
		checkbox.setChecked(layer.isDisplayed());
		checkbox.addEventListener("onCheck", new VisibilityToggleEventListener());
		checkbox.setParent(listcell);
		checkbox.setTooltiptext("Hide");
		
		Label label = new Label(LayerUtilities.chompLayerName(layer.getName()));
		label.setParent(listcell);
		listcell.setParent(item);

		// dnd list reordering support
		item.addEventListener("onDrop", new ActiveLayerDNDEventListener());
		item.setDraggable("true");
		item.setDroppable("true");
		
		// bind to the ActiveLayer instance (we readback later)
		item.setValue(layer);
	
		// simple description for tooltip
		label.setTooltiptext(LayerUtilities.getTooltip(layer.getName(),layer.getDescription()));
		
		label.setStyle("float:left;");
		checkbox.setStyle("float:left;");
		

		
		/* show the legend graphic when the user hovers over the pallette icon
		 * if 
		 */
		Image remove = new Image(Config.getLang("layer_remove_icon"));
		remove.addEventListener("onClick", new ActiveLayersRemoveEventListener());
		remove.setParent(listcell);
		remove.setStyle("float:right;");
		remove.setTooltiptext("remove layer");
		
		if (layer.isDefaultStyleLegendUriSet()) {
			/* animated layers get a "special" key icon
			 * Since setting this icon will only happen if style 
			 * legends are available, it is an assumption that
			 * only layers with legends will be animatable, so 
			 * if you have a layer with no legend that doesn't 
			 * display the special animation key icon, that's 
			 * why
			 */
			Image legend;
			if (layer.isCurrentlyAnimated()) {
				// layer is currently being animated (implies supports animation)
				legend = new Image(Config.getLang("map_legend_animated_icon"));
			}
			else if (layer.isSupportsAnimation()) {
				// layer is supports animation but is not currently animated
				legend = new Image(Config.getLang("map_legend_animatable_icon"));
			}
			else {
				// just a plain layer
				legend = new Image(Config.getLang("map_legend_icon"));
			}
			
			/* hack to get things to align properly - want everything
			 * floating left except the image which floats right 
			 */

			legend.setStyle("float:right;");
			
			legend.setParent(listcell);
			
			// hover a tooltip image over the icon
			Popup popup = (Popup) Executions.createComponents("/WEB-INF/zul/LegendPopup.zul", legend.getRoot(), null);	
			popup.addEventListener("onOpen", new LegendTooltipOpenEventListener(layer));
			legend.setTooltip(popup);
			legend.addEventListener("onClick", new LegendClickEventListener(layer));

		}		
	}
}
