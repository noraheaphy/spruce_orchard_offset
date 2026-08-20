# Spruce orchard genomic offset paper
Scripts and data to reproduce analyses in spruce orchard genomic offset paper.

**Data**
- `clim_sprucePops.csv` -- 11 climate variables (1961-1990) for 64 sampled red spruce populations
- `snpScores_sprucePops.csv` -- allele frequencies of 240 candidate SNPs for 64 sampled red spruce populations
- `currClim_sprucePres_1960_1990.csv` -- 11 climate variables (1961-1990) for all geographic grid cells (2.5 arc-min) with red spruce presences (8821 lat/lons)
- `futClim_2071_SSP245.csv.zip` -- 11 climate variables (2071-2100, SSP2-45) for all geographic grid cells (2.5 arc-min) within study area
- `US_states.shp` -- shapefiles for US states for clipping presences to study area

**Scripts**
- `01_rawGenomicOffsets.R` —- uses gradient forest modeling to calculate raw genomic offsets, both spatial (between all current red spruce presences) and spatiotemporal (between current donors and future recipient sites)
- `02_standardizedOffsets.R` —- standardizes spatiotemporal offsets based on empirical distribution of current spatial offsets, applies threshold to binarize offsets into acceptable and unacceptable transfers, calculates donor importance and recipient importance
- `03_recImp_donorsScions.R` -- calculates recipient importance for top candidate donors as well as actual scion collected for orchard
- `04_alternateSelectionMethods.R` -- compares recipient importance based on genomic offset-based candidate donor selection vs. alternative methods of seed sourcing (random, local to orchard site, within same seed zone, one seed zone north, from within seed zones containing top genomic offset candidates)
