bed_file = open('D://Dropbox/bioinfo/final/hg19_bed.txt', 'r')
bed_line = (x for x in bed_file.read().splitlines())
bed_file.close()
bed_element = (line.strip().split('\t') for line in bed_line)
exon_dic = {'chromosome':[], 'start':[], 'end':[]}
intron_dic = {'chromosome':[], 'start':[], 'end':[]}
for i in bed_element:
    gene_st = int(i[1])
    gene_ed = int(i[2])
    exon_num = int(i[9])
    exon_size = [int(x) for x in i[10].split(',')[:-1]]
    exon_start = [int(x) for x in i[11].split(',')[:-1]]
    for j in range(exon_num):
        if (j < (exon_num-1)):
            exon_dic['chromosome'].append(i[0])
            exon_dic['start'].append(gene_st + exon_start[j])
            exon_dic['end'].append(gene_st + exon_start[j] + exon_size[j])
            intron_dic['chromosome'].append(i[0])
            intron_dic['start'].append(gene_st + exon_start[j] + exon_size[j] +1)
            intron_dic['end'].append(gene_st + exon_start[j+1] -1)
        else:
            exon_dic['chromosome'].append(i[0])
            exon_dic['start'].append(gene_st + exon_start[j])
            exon_dic['end'].append(gene_st + exon_start[j] + exon_size[j])

from pandas import DataFrame
exon_df = DataFrame(exon_dic, columns=['chromosome', 'start', 'end'])
intron_df = DataFrame(intron_dic, columns=['chromosome', 'start', 'end'])

ndc_chr = ['chr'+str(x) for x in range(1,23)+['X','Y']]
exon_df = exon_df[exon_df['chromosome'].isin(ndc_chr)]
intron_df = intron_df[intron_df['chromosome'].isin(ndc_chr)]

exon_df.to_csv('exon.csv', index=False)
intron_df.to_csv('intron.csv', index=False)