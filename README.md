# Wz_CM
This R code can be used to extract daily vertical current velocity as a weighted mean of a user specified maximum depth to marine occurence data (e.g., satellite tracking data or sightings data), using Copernicus Marine data (https://data.marine.copernicus.eu/products) already downloaded to a hard drive using the Copernicus Marine Toolbox (https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-toolbox-introduction). The code presented here is a modified version of code available in the "Remora" R package available on Github (https://github.com/IMOS-AnimalTracking/remora) which will use the Bluelink BRAN2020 product (https://research.csiro.au/bluelink/bran2020-data-released/) as data source. However, BRAN2020 is only available until the end of 2023. The code presented here can be used to extract the latest 2024 data or older data, depending on the availability of COpernius Marine data. The code here was written to work for data structure of the Global Ocean Physics Analysis and Forecast product with product ID: GLOBAL_ANALYSISFORECAST_PHY_001_024 (https://doi.org/10.48670/moi-00016), and the spcific Datset ID: cmems_mod_glo_phy-wcur_anfc_0.083deg_P1D-m.
