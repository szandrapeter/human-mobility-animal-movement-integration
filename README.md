## **Jupyter Notebooks and R scripts used to integrate human mobility and animal movement data**

## **Integrating human mobility and animal movement data reveals complex space-use between humans and white-tailed deer in urban environments**

### **Authors**

Szandra A. Péter¹, Travis Gallo², Jennifer Mullinax², Amira Roess³, Gabriela Palomo-Munoz²˒⁴, Taylor Anderson¹
1.	Department of Geography and Geoinformation Science, George Mason University, Fairfax, VA, 22030, USA
2. 	Department of Environmental Science and Technology, University of Maryland, College Park, MD, 20742, USA
3.	Department of Global and Community Health, George Mason University, Fairfax, VA, 22030, USA 
4.	Biology Department, Universidad del Valle de Guatemala, Guatemala City, Guatemala (current affiliation)

---

### **File Descriptions**

#### **Jupyter Notebooks**

- `Scripts/movement_buildings_poptotal.ipynb` 
  Computes and processes movement metrics (visit counts, popularity by hour) for deer and humans, and calculates the residential and commercial building area, as well as the estimated nighttime population, across a hexagonal spatial tessellation. Also merges park hexagons and aggregates the corresponding values. Produces `visit_counts_201801_202001.shp` which is read into `2025-03-24_data_processing_hex.R`.
  
- `Scripts/GEE_temperature_download.ipynb` 
  Downloades the publicly available daytime MODIS Land Surface Temperature (LST) data using Google Earth Engine. 

- `Scripts/figure3.ipynb` 
  Computes mean popularity by hour of deer and humans across different landscape types. The resulting data were used to produce Figure 3.

- `Scripts/figure4.ipynb` 
  Computes seasonal mean popularity by hour of deer and humans across different landscape types. The resulting data were used to produce Figure 4.

#### **R Scripts**

- `Scripts/package_load.R` 
  Loads needed libraries. Source this script at the beginning to ensure all necessary packages and functions are available.

- `Scripts/2025-03-24_data_processing_hex.R` 
  Processes landscape-level variables for the hexagons. Produces `2025-03-24_full_data.rds` which is read into `2025-03-24_human_deer_combined_nb_lasso_hex.R` and `2025-03-24_nb_offset_deer_nimble_hex.R`. Note that this script requires several spatial datasets that are publicly available, but we do not provide in this repository. Sources are indicated in the script and/or the manuscript.

- `Scripts/2025-03-24_human_deer_combined_nb_lasso_hex.R` 
  Independent negative binomial models for deer and humans. Produces `2025-03-24_CombinedAnalysis_DeerModelResults.rds` and `2025-03-24_CombinedAnalysis_HumanModelResults.rds`.
   
- `Scripts/2025-03-24_nb_offset_deer_nimble_hex.R` 
  Deer-only negative binomial model using the same independent variables as `2025-03-24_human_deer_combined_nb_lasso_hex.R`, with human commercial activity included as an additional variable. Produces `2025-03-24_model_results_nb_lasso.rds`.

---

#### **Data Files**

- `Data/deer_GPS_2018_2020.zip`  
  Contains the cleaned shapefile of hourly GPS-collar data from deer in Howard County, Maryland. This is a subset of data from Roden-Reynolds *et al.* (2022). Each record corresponds to an individual-deer presence at a given hour and location. Read into `movement_buildings_poptotal.ipynb`.

| Field Name         | Data Type | Description                                             |
|--------------------|-----------|---------------------------------------------------------|
| OBJECTID           | Object ID | Unique identifier for each record                       |
| Shape              | Geometry  | Spatial geometry of the grid cell                       |
| event_id           | Double    | Event identifier of the GPS recording                   |
| loc_long           | Double    | Longitude coordinate                                    |
| loc_lat            | Double    | Latitude coordinate                                     |
| tag_loc_id         | Text      | Tag local identifier of the deer                        |
| ind_loc_id         | Long      | Unique, individual local identifier of the deer         |
| gmt_dt             | Text      | Timestamp in UTC                                        |

- `Data/visit_counts_201801_202001/visit_counts_201801_202001.shp`  
  Contains aggregated visit counts of deer and humans for grid cells after merging park hexagons. Each record includes the monthly visit counts of deer and humans within a specific grid cell, which may be a single hexagon or a merged area of multiple hexagons. For each record, the POI count, residential and commercial building area, estimated total nighttime population, and the number of original hexagons merged (>1 if park hexagons were combined) are also included. Read into `2025-03-24_data_processing_hex.R`.

| Field Name     | Data Type | Description                                                                     |
|----------------|-----------|---------------------------------------------------------------------------------|
| OBJECTID       | Object ID | Unique identifier for each record                                               |
| Shape          | Geometry  | Spatial geometry of the grid cell                                               |
| GRID_ID_p      | Text      | Unique grid cell identifier                                                     |
| Shape_Length   | Double    | Perimeter length of the grid cell                                               |
| Shape_Area     | Double    | Area of the grid cell                                                           |
| Jan18_deer - Jan20_deer | Double  | Monthly deer visit counts (e.g., `Jan18_deer` = deer visits in Jan 2018) |
| Jan18_hum - Jan20_hum | Double  | Monthly human visit counts (e.g., `Jan18_hum` = human visits in Jan 2018)  |
| POI_count      | Double    | Number of points of interest (POIs) in the grid cell                            |
| res_m2         | Double    | Residential building area (m²) within the grid cell                             |
| com_m2         | Double    | Commercial building area (m²) within the grid cell                              |
| pop_total      | Double    | Estimated nighttime population within the grid cell                             |
| GRID_count     | Long      | Total number of original hexagons merged                                        |

- `Data/hexagonal_tessellation/hexagonal_tessellation.shp`  
  Spatial data of the hexagonal tessellation used in the study. Read into `movement_buildings_poptotal.ipynb`.

- `Data/2025-03-24_full_data.rds`  
  Processed data used in the negative binomial models. Read into `2025-03-24_human_deer_combined_nb_lasso_hex.R` and `2025-03-24_nb_offset_deer_nimble_hex.R`.
  
---

### **Notes**
- The Jupyter Notebooks (except for `GEE_temperature_download.ipynb`) were run within ArcGIS Pro and include `arcpy` functions. These functions are comparable to those available in other Python packages and could be adapted accordingly.

---

### **Reference**
Roden-Reynolds, P., Kent, C. M., Li, A. Y., & Mullinax, J. M. Patterns of white-tailed deer movements in suburban Maryland: implications for zoonotic disease mitigation. *Urban Ecosyst.* **25**, 1925–1938 (2022).