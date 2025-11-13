---
title: "Notebook 1"
layout: single
sidebar:
  nav: "docs"
---


# Notebook 1: Air Pollution Data Preparation

## README

### Overview
This notebook conducts the data processing and comprehensive analysis of nitrogen dioxide (NO₂) pollution using [Sentinel-5P data](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S5P_NRTI_L3_NO2), with a focus on **Ethiopia (Addis Ababa)** and **Iraq (Baghdad)**. It covers the full workflow, including data retrieval, preprocessing, aggregation, and visualisation.

In addition to NO₂, other key pollutants - **sulfur dioxide (SO₂)**, **carbon monoxide (CO)**, and **ozone (O₃)** - are also processed and briefly explored for inter-pollutant relationships.

### Objective
- **Clean & aggregate** raw Sentinel-5P NO₂ to mesh grids.  
- **Visualize** spatial–temporal patterns (static maps & GIFs).  
- **Analyse** pollutant dynamics with  
  1. **Time-series plots** (daily & monthly trends), and  
  2. **Correlation matrix** across NO₂, SO₂, CO, and O₃.  

### Workflow
1. **Data Retrieval**  
   - Download NO₂, SO₂, CO, O₃ via Earth Engine Python API.

2. **Data Processing**  
   - Gap-fill (e.g., cloud masking) and clip to city boundaries.  
   - Aggregate NO₂ to mesh grids; optionally resample to monthly means.

3. **Exploratory Analysis**  
   - **Time-Series Analysis**: plot daily & monthly pollutant trends.  
   - **Correlation Analysis**: compute Pearson/Spearman matrix and visualise as heat-map.

4. **Visualisation**  
   - Compare raw vs. filled rasters.  
   - Overlay NO₂ on urban grids.  
   - Export animated GIFs for NO₂ evolution.

### Data Process Pipeline

This notebook processes the air pollution data downloaded in *appendix_preparation.ipynb* through the following steps:

- **(1) Filling Missing Value**: Spot the missing values in raster and replenish them using iterative filling, using **mean** of the neighbour raster as the replenish value.

- **(2) Clipping to Region**: Clipping the data to the interested area, and output the filled raster.

- **(3) Aggregation**: Import the generated mesh and aggregate the raster to the mesh level.

Step 2 and 3 are realised by selecting and aggregating the data within the mesh grid. 

### Outputs
| Type | File | Description |
|------|------|-------------|
| Raster | `no2_filled_rasters.tif` | Gap-filled NO₂ rasters |
| Vector/Tabular | `no2_mesh_data.gpkg`| Grid-level NO₂ statistics |
| Figures | `visual_comparison.png` | Raw vs. filled snapshot |
| Figures | `no2_timelapse.gif` | Time-lapse of NO₂ |
| Figures | `pollutant_timeseries.png` | Multi-pollutant time-series plots |
| Figures | `pollution_corr_heatmap.png` | Correlation matrix heat-map |

*Note: no2 above can be changed to other types of air pollution*

### Notes
- The NO₂ filling step is computationally heavy and may take 8+ hours per city for full-year data.

- Data quality is significantly improved after gap-filling, especially in cloud-obscured areas.