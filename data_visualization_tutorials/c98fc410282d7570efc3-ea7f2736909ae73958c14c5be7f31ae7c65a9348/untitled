#!/usr/bin/env python

"""
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.

1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::

    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg

2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.

"""

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """

    Yields BrokenBarHCollection of features that can be added to an Axes
    object.

    Parameters
    ----------

    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.

    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection

    height : float
        Height of each BrokenBarHCollection

    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print chrom
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


# Height of each ideogram
chrom_height = 1 #the vertical length of the chromosome bar

# Spacing between consecutive ideograms
chrom_spacing = 1 #spacing between consecutive chromosome entries

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.4 #the vertical length of the blue lines (gene track)

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.1 

# Width, height (in inches)
figsize = (6, 8) #size of the window in inches

# Decide which chromosomes to use
chromosome_list = ['chr%s' % i for i in range(1, 23) + ['M', 'X', 'Y']] #

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0 #represents the bottom of the chrom bar
chrom_ybase = {} #represents mapping of  chrom --> y coord of bottom of chrom bar
gene_ybase = {}  #represents mapping of chrom --> y coord of bottom of blue lines
chrom_centers = {} #represents mapping of chrom --> y coord of center of chrom bar

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]: #::-1 means reverse
    chrom_ybase[chrom] = ybase #sets the y coord of the chrom bar bottom to the value of ybase
    chrom_centers[chrom] = ybase + chrom_height / 2. #sets the y coord of the center of the bar to halfway point
    gene_ybase[chrom] = ybase - gene_height - gene_padding #sets the y coord of bottom of blue lines to 0.5 below chrom bar
    ybase += chrom_height + chrom_spacing #increments ybase for the next chromosome by 2

# Read in cytoBandIdeo.txt, downloaded from UCSC Table Browser
# cytoBandIdeo.txt is a file that contains the chromosome, the start and end of each chromosome, the name, and the color of the stain
# ideo is a DataFrame
ideo = pandas.read_table( 
    'cytoBandIdeo.txt', #the file we are reading in
    skiprows=0,#skip the first row
    names=['chrom', 'start', 'end', 'name', 'gieStain']#names of the columns
)

# Filter out chromosomes not in our list
ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)] #ask Denise

# Add a new column for width
ideo['width'] = ideo.end - ideo.start #subtracts end value - start value and puts into width col

# Colors for different chromosome stains
color_lookup = { # A dictionary that maps from color string --> rgb values
    'gneg': (1., 1., 1.),
    'gpos25': (.6, .6, .6),

    'gpos50': (.4, .4, .4),
    'gpos75': (.2, .2, .2),
    'gpos100': (0., 0., 0.),
    'acen': (.8, .4, .4),
    'gvar': (.8, .8, .8),
    'stalk': (.9, .9, .9),
}

# Add a new column for colors
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x]) #=
print ideo['colors'] #remove later

# Same thing for genes
genes = pandas.read_table(
    'ucsc_genes.txt',
    names=['chrom', 'start', 'end', 'transcript', 'name'],
    usecols=range(5))
genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
genes['width'] = genes.end - genes.start
genes['colors'] = '#2243a8'


fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
    ax.add_collection(collection)

# ...and the gene data
print("adding genes...")
for collection in chromosome_collections(
    genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
):
    ax.add_collection(collection)

# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.axis('tight')
plt.show()