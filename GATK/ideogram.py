import matplotlib.pyplot as plt
import pandas

#chrom changes file
col_names = ['Read ID', 'chromosome', 'Y Starting Base Position', 'Percent Match', 'Second Position', 'ChangeTo', 'Starting Base Position', 'Percent Match2', 'Second Position2']
chrom_changes = pandas.read_table( 
    'LP6005636-DNA_H02.chrom_changes.txt', #the file we are reading in
    names= col_names)#names of the columns
 
chrom_changes.set_index("Read ID", drop=True, inplace = True) #sets index to first col
chrom_changes.head()
#chrom_changes.loc[chrom_changes['Second Position'] != "="] ask denise about this

#no match file
col_names2 = ['Read ID', 'chromosome', ' Y Starting Base Position', 'Percent Match', 'Second position', 'Genome w/out y' ]
no_match_reads = pandas.read_table( 
    'LP6005636-DNA_H02.no_match.sorted.txt', #the file we are reading in
    names= col_names2)#names of the columns
 
no_match_reads.set_index("Read ID", drop=False, inplace = True) #sets index to first col
no_match_reads.head()
#no_match_reads.loc[no_match_reads['Second position'] != "="] checking to make sure there aren't extra dots to plot

# removing duplicate reads on the y from no match file
y_locations = no_match_reads[' Y Starting Base Position']
unique_y_locations = set(y_locations)
#for index, position in y_locations.iteritems():
    #if position not in unique_y_locations:
        #unique_y_locations.append(position)


 #plotting the y chromosome ideogram
ax1 = plt.subplot(211)
ax1.set_title("An ideogram of read locations for chromsome Y ~ no matches file")
ax1.set_ylabel("Chromosome")
ax1.set_xlabel("Starting BP location of read")
ax1.set_yticklabels(["", "", "Chr Y", "", ""])
plt.axis([0.25 * 10**7, 0.7 * 10**7, 0, 2])

for position in unique_y_locations:
    ax1.plot(position, 1, 'b.')
#plt.show()


# lengths of chromosomes taken from https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
# lengths listed in order from 1-22, X, M, Y
chr_lengths = [249250621, 243199373,198022430,191154276, 180915260, 
               171115067, 159138663, 146364022, 141213431, 
               135534747, 135006516, 133851895, 115169878, 107349540, 
               102531392, 90354753, 81195210, 78077248, 59128983,
               63025520, 48129895, 51304566, 155270560, 16571, 59373566]  


 #creating a dictionary that links chrom labels to y coordinates
chrom_list = ['chr%s' % i for i in range(1, 23) + ['X', 'M', 'Y']]
chrom_dict = {}
for i, chrom in enumerate(chrom_list):
    chrom_dict[chrom] = i + 1


 #plotting the lines for the chrom changes ideogram
ax2 = plt.subplot(212)
ax2.set_title('Ideogram for change chroms file')
ax2.set_xlabel('Starting Base Position')
ax2.set_ylabel('Chromosome')
ax2.set_yticklabels(chrom_list)
for i in range(0, 25):
    ax2.hlines(i + 1, 0, chr_lengths[i], linestyles='solid')
#plt.show()


#iterating over each row of chrom changes so that points can be plotted onto ideogram
for index, row in chrom_changes.iterrows():
    #plotting chrY dot
    ax2.plot(row['Y Starting Base Position'], chrom_dict['chrY'], 'b.')
    #plotting changed read dot
    ax2.plot(row['Starting Base Position'], chrom_dict[row['ChangeTo']], 'b.')
plt.show()   