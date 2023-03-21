import csv

#  Receives a .fasta file of gene CDS's and filters out any CDS that does not begin with ATG or end with a stop codon
#  The output file is used as input for third party software CodonW which computes GC3 content for each gene

csv_data = []  # Filtered list of CDS output
CDS_rows = []  # List of the row number that each CDs begins at in the input file

with open('Rhoto_IFO0880_4_GeneCatalog_CDS_20170509.fasta', 'r') as Source:
    report_reader = csv.reader(Source, delimiter='\t')
    report_reader_list = list(report_reader)
    for idx, row in enumerate(report_reader_list):
        try:
            if row[0][0] == '>':
                CDS_rows.append(idx)  # Identifies & saves locations of title rows for each gene in .fasta file to
                # help find the transcriptID, and start and stop codons
        except IndexError:
            continue
    for idx, locus in enumerate(CDS_rows):
        try:
            if report_reader_list[locus + 1][0][0:3] == 'ATG' and (
                    report_reader_list[CDS_rows[idx + 1] - 1][0].endswith('TGA') or
                    report_reader_list[CDS_rows[idx + 1] - 1][0].endswith('TAG') or
                    report_reader_list[CDS_rows[idx + 1] - 1][0].endswith('TAA')):
                # Checks if each gene starts & ends with start & stop codon, if so adds gene to filtered list
                if report_reader_list[locus][0][25] == '|':  # The relevant portion of the gene title is either 4 or 5
                    #  digits followed by the | character
                    csv_data.append('>' + report_reader_list[locus][0][21:25])
                else:
                    csv_data.append('>' + report_reader_list[locus][0][21:26])
                for x in range(locus + 1, CDS_rows[idx + 1]):
                    csv_data.append(report_reader_list[x][0])
        except IndexError:
            continue

output = open('CDS_processed.fasta', 'w')
for row in csv_data:
    output.write(row)
    output.write('\n')
output.close()
