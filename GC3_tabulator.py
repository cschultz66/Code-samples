import csv

#  Receive CodonW output with transcriptID-GC3-content pairs and use .gff3 gene annotation file to generate list of
#  each gene in a given assembly scaffold along with that gene's start and end locus and GC3 content

gene_loci = []  # Output list of genes by transcriptID, start locus, end locus, GC3 content
GC3_parsed = []  # List of genes by transcriptID, GC3 content after processing to remove spaces and tabs from CodonW
# output

with open('Rhoto_IFO0880_4.FrozenGeneCatalog.gff3', 'r') as Source:
    report_reader = csv.reader(Source, delimiter='\t')
    annotationsList = list(report_reader)

with open('CDS_processed.out', 'r') as Source:
    GC3_reader = csv.reader(Source, delimiter='\t')
    GC3_reader_list = list(GC3_reader)

for row in GC3_reader_list:  # Remove spaces and tabs from first column of CodonW output
    GC3_parsed.append([row[0][0:row[0].find(' ')], row[1]])

for row in annotationsList:  # Checking each object in the annotation file
    try:
        if row[0] == 'scaffold_18' and row[2] == 'gene':  # If the object is a gene on the chosen scaffold
            GC3content = 'NA'
            for GC3row in GC3_parsed:  # Find the gene with matching transcriptID from CodonW
                if GC3row[0] == row[8][row[8].find('transcriptId=') + 13:len(row[8])]:
                    GC3content = GC3row[1]
            if GC3content != 'NA':  # If a match was found, save the transcriptID, start locus, end locus, GC3 content
                gene_loci.append([row[8][row[8].find('transcriptId=') + 13:len(row[8])], row[3], row[4], GC3content])
    except IndexError:
        continue

with open('scaffold_18_GC3_v2.txt', 'w') as output:
    report_writer = csv.writer(output, delimiter='\t')
    for row in gene_loci:
        if row[0] != '':
            report_writer.writerow(row)
