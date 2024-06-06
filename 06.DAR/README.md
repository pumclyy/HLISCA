# Analysis of smoking-responsive cCREs

We employed the "FindMarkers" function in Seurat to identify peaks of differential accessibility between the smoking and non-smoking states for each cell type.

Initially, we divided each 23 cell type into two clusters: smoker cluster and non-smoker cluster. Subsequently, we utilized the MACS2 algorithm to identify significant peaks within the 46 resultant clusters.

We performed findmarkers between smoker cluster and non smoker cluster of each cell type to identify smoking-responsive cCREs. Then we annotated smoking-responsive cCREs using ChIPseeker to analysis these cCREs.
