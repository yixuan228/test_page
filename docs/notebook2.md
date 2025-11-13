---
title: "Notebook 2"
layout: single
sidebar:
  nav: "docs"
---


# Notebook 2: Related Features Preparation

## README

### Overview
This notebook explores key spatial datasets related to human activities and infrastructure that may contribute to or correlate with nitrogen dioxide (NO₂) pollution in two case study cities: Addis Ababa, Ethiopia and Baghdad, Iraq. Including:

   - **Point of Interest and Road Network Data**, from Open Street Map for [Ethiopia](https://download.geofabrik.de/africa/ethiopia-latest-free.shp.zip) and [Iraq](https://download.geofabrik.de/asia/iraq-latest-free.shp.zip).
   - **Population Data**: [WorldPop](https://data.humdata.org/search?q=Spatial+Distribution+of+Population+&ext_search_source=main-nav) for population density.
   - **Night Time Light Data**, from [NASA Black Marble](https://blackmarble.gsfc.nasa.gov) for night time light.
   - **Cloud Cover Data**, from [MODIS/Aqua Surface Reflectance](https://lpdaac.usgs.gov/products/myd09a1v061/) for cloud cover.
   - **Land Surface Temperature Data**, from [MODIS/Terra Land Surface Temperature](https://lpdaac.usgs.gov/products/mod11a1v061/) for land surface temperature.
   - **Traffic Congestion Data**, from WorldBank internal source, **only available in Baghdad**.

   - Auxiliary dataset: **Boundary data**, from [HDX](https://data.humdata.org/) for [Ethiopia](https://data.humdata.org/dataset/cod-ab-eth) and [Iraq](https://data.humdata.org/dataset/cod-ab-irq).


### Objective
The main objective of this notebook is to:

- Finish the data processing of different datasets, including missing data filling, clipping to interested area, aggregation and basic visualization.

- Provide mesh-based visualizations of urban-scale features such as population density, road length, and POI distribution, to support future modeling

### Workflow

1. **Data Acquisition**
   - Downloads datasets from multiple sources listed above.

2. **Data Processing**
   - Clips all spatial layers to city-level boundaries. 
   - Aggregates various indicators into uniform mesh grids for spatial comparison.

3. **Visualisation**
   - Generates static maps of population density and road infrastructure distribution.
   - Visualises spatial variation in contributing factors to support exploratory analysis.

4. **Export Results**
   - Save enriched mesh grids with new attributes for modeling and further analysis.

### Data Process Pipeline

The datasets in this notebook are processed through the following steps(similar to the process in notebook1):

- **(1) Filling Missing Value**: Spot the missing values in raster and replenish them using iterative checking, using **mean** of the neighbour raster as the replenish value.

- **(2) Clipping to Region**: Clipping the data to the interested area, and output the filled raster.

- **(3) Aggregation**: Import the generated mesh and aggregate the raster to the mesh level.

- **(4) Visualization**: Visualize the result as static figure [*.png*] or dynamic figure [*.gif*].

### Outputs

- **Processed GeoPackages**: City-level road, POI, and infrastructure layers.
- **Mesh-Based Aggregated Tables**: Feature counts and densities per grid cell.
- **Static Visualisations**: Color maps for all contributing factors to aid spatial interpretation.
- **Harmonised Datasets**: All outputs aligned spatially and ready for merging with NO₂ concentration data.