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

   
