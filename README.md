# Palm_Oil_Recovery



## Study info:

There are 16 samples from each of 5 plots. All plots are from a palm oil plantation that is undergoing [rewilding](https://alittlewild.com/#rec523988346).

**Plots:**

  - P5 is oil palm, still farmed
  - P4 is an area where all the oil palm trees have been removed – has been like this for 1yr (left fallow)
  - P3 a managed plot that has been managed/plated for three years (e.g., three years since clearing)
  - P2 a managed plot that has been managed/plated for two years (e.g., two years since clearing)
  - P1 a managed plot that has been managed/plated for one years (e.g., one years since clearing)

Using Plot (P3) three as an example, it was palm oil originally, was cleared, and has been managed for three years.

[**This link**](https://www.google.com/maps/d/viewer?mid=1Jje6npngtSd5nsdgyZiSPmr96sGx1B8&ll=1.6815873025227588%2C103.83737993658451&z=20) shows where the plots are and has a few pics of each plot. 


**In all plots they follow the same regime:**

  - Palm oil production
  - Area cleared and left fallow for a year
  - Then they plant X (will find our exactly what) in year 1
  - Then Y in year 2
  - And Z in year 3

All started off as palm oil, then spent a year clear and left alone before planting X, Y then Z 

___

## Data info:

Raw data is on the SRA: [PRJNA1117193](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1117193)

The data is full length 16S and both regions of fungal ITS

Taxonomy assigned against the following databases:

  - Bacteria -- [SILVA v 138.1](10.5281/zenodo.4587946)
  - Eukaryotes -- [Eukaryome v 1.9](https://eukaryome.org/)

___


**Starting hypotheses:**

  - 1. Microbial diversity increases with rewilding age and we see corresponding shifts in community structure that influences the chemical properties of the soil. Or something along those lines.
  - 2. Network complexity in bacteria and fungi increases with time


**Initial analyses:**

  - 1, Standard Plots of diversity for bacteria and fungi
  - 2, NMDS plots Bacteria and fungi – I would be surprised if they do not cluster by plot.
  - 3, CornCob and all the other methods we use, differential abundance for bacteria and fungi – what is driving
  - 4, Hopefully with these longer read length we can assign some functions to those that are differentially abundant and these functions relate to the measured pH, C, N etc.
  - 5, Heatmaps of top 20 taxa for both, maybe we can resolve to genus or even species now?!!
  - 6, Co-occurrence networks, maybe even a bit of complexity work – I suspect the palm oil is the least complex, whereas the plot that has been managed for three years is the most complex
  - 7, Keystone species in the networks like in this paper (fig 3, https://doi.org/10.1007/s00248-024-02372-5)
  - 8, Does soil chemical property influence community structure?
  - 9, We might be able to get at community assembly a bit like fig 2 in https://doi.org/10.1007/s00248-024-02372-5, instead of tree species we use plot
  - 10, anything else you want to try? Any ideas from the restoration work you’ve done in the US etc?

 
Convert "assembly over time taxa plots" to phylogeny: 4 trees, left to right; highlight branches that appear in that year or previous years; expanding tree of diversity!
 
 
 
 