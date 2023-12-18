# Spatial Temporal Distribution of b-value (STDB)
A versatile b-value calculator. Used for calculating time and space variations of b-values, it includes eight Mc calculation methods and four b-value calculation methods.
Additional functions include synthetic seismic catalogs for seismic events, and more.

STDB is an extremely simple b-value calculator that includes most algorithms on the market, and the grouping between different algorithms is flexible and efficient.
Just import a seismic catalog.

# Earthquake Catalog Format (the following numbers represent column numbers):
1.Year
2.Month
3.Day
4.Hour
5.Minute
6.Latitude
7.Longitude
8.Magnitude

# When estimating the b-value using the STDB calculator, refer to:
Weicheng Gong, Huayuan Cheng, Yajing Gao, Qing Li, and Yunqiang Sun (2023). Spatial-Temporal Distribution of b-value in the Eastern Tibetan Plateau. (in revision)

# STDB calculator's main features:
1. Calculate the b-value for the entire catalog
2. Compute the b-value as a function of time
3. Determine the b-value as a function of space
4. Generate a synthetic earthquake catalog (additional feature)
5. Perform iterative sampling on the catalog (additional feature)

# STDB calculator involves the following Mc calculation methods:
1. The Maximum Curvature technique (MAXC) (Wiemer & Wyss, 2000)
2. The goodness-of-fit test (GFT) (Wiemer & Wyss, 2000)
3. The Mc by b-value stability (MBS) (Cao & Gao, 2002)
4. The improved Mc by b-value stability (MBS-WW) (Woessner & Wiemer, 2005)
5. The entire-magnitude-range (EMR) (Woessner & Wiemer, 2005)
6. The median-based analysis of the segment slope (MBASS) (Amorese, 2007)
7. The harmonic mean of the magnitudes (MAH) (Godano, 2017)
8. The Normalized Distance test (NDT) (Lombardi, 2021)

# STDB calculator involves the following b-value calculation methods:
1. linear least squares regression (LSR) (Pacheco et al., 1992)
2. maximum likelihood estimation (MLE) (Aki, 1965)
3. B-Positive (van der Elst, 2021)
4. K-M slope (KMS) (Li et al., 2023; Telesca et al., 2013)

# Functions involved in the STDB calculator's synthetic seismic catalog:
1. OK1993 model Ogata & Katsura (1993)
2. WW2005 model Woessner & Wiemer (2005)

# Other references
Aki, K. (1965). Maximum likelihood estimate of b in the formula logN=a-bM and its confidence limits. Bulletin of the Earthquake Research Institute, University of Tokyo, 43, 237-239.
Amorese, D. (2007). Applying a change-point detection method on frequency-magnitude distributions. Bulletin of the Seismological Society of America, 97(5), 1742-1749. https://doi.org/10.1785/0120060181
Cao, A., & Gao, S. S. (2002). Temporal variation of seismic b-values beneath northeastern Japan island arc. Geophysical Research Letters, 29(9), 1334. https://doi.org/10.1029/2001GL013775
Godano, C. (2017). A new method for the estimation of the completeness magnitude. Physics of the Earth and Planetary Interiors, 263, 7-11. https://doi.org/10.1016/j.pepi.2016.12.003
Gulia, L., & Wiemer, S. (2019). Real-time discrimination of earthquake foreshocks and aftershocks. Nature, 574(7777), 193-199. https://doi.org/10.1038/s41586-019-1606-4
Li, L., Luo, G., & Liu, M. (2023). The K-M Slope: A Potential Supplement for b-Value. Seismological Research Letters. https://doi.org/10.1785/0220220268
Lombardi, A. M. (2021). A Normalized Distance Test for Co‐Determining the Completeness Magnitude and b‐Value of Earthquake Catalogs. Journal of Geophysical Research: Solid Earth, 126(3). https://doi.org/10.1029/2020jb021242
Ogata, Y., & Katsura, K. (1993). Analysis of temporal and spatial heterogeneity of magnitude frequency distribution inferred from earthquake catalogues. Geophysical Journal International, 113(3), 727-738. https://doi.org/10.1111/j.1365-246X.1993.tb04663.x
Pacheco, J. F., Scholz, C. H., & Sykes, L. R. (1992). Changes in frequency-size relationship from small to large earthquakes. Nature, 355(6355), 71-73. https://doi.org/10.1038/355071a0
Telesca, L., Lovallo, M., Ramirez-Rojas, A., & Flores-Marquez, L. (2013). Investigating the time dynamics of seismicity by using the visibility graph approach: Application to seismicity of Mexican subduction zone. Physica A: Statistical Mechanics and its Applications, 392(24), 6571-6577. https://doi.org/10.1016/j.physa.2013.08.078
van der Elst, N. J. (2021). B‐Positive: A Robust Estimator of Aftershock Magnitude Distribution in Transiently Incomplete Catalogs. Journal of Geophysical Research: Solid Earth, 126(2). https://doi.org/10.1029/2020jb021027
Wiemer, S., & Wyss, M. (2000). Minimum magnitude of completeness in earthquake catalogs: Examples from Alaska, the western United States, and Japan. Bulletin of the Seismological Society of America, 90(4), 859-869. https://doi.org/10.1785/0119990114
Woessner, J., & Wiemer, S. (2005). Assessing the quality of earthquake catalogues: Estimating the magnitude of completeness and its uncertainty. Bulletin of the Seismological Society of America, 95(2), 684-698. https://doi.org/10.1785/0120040007

   
