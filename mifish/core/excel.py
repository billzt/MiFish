import xlsxwriter
import json
import os

def read_stat(workdir:str, sample2type:dict):
    workbook = xlsxwriter.Workbook(f'{workdir}/MiFishResult/read_stat.xlsx')
    worksheet = workbook.add_worksheet()
    header_format = workbook.add_format({'bg_color': 'yellow', 'font_name': 'Arial', 'align':'center'})
    headers = 'Sample Name	Read Type	No. Reads	No. Reads after Quality Filter	No. Assembled Reads	No. Reads Without N	No. Reads Fit Length	No. Unique Reads	No. Denoised Haploids	No. Species'.split('\t')
    worksheet.set_column(0,0,30)
    worksheet.set_column(1,1,15)
    worksheet.set_column(2,len(headers)-1,20)
    for col in range(len(headers)):
        worksheet.write(0, col, headers[col], header_format)
    
    body_format = workbook.add_format({'font_name': 'Arial'})
    row = 1
    for (sample_name, input_type) in sample2type.items():
        stat_file = f'{workdir}/MiFishResult/Sample-{sample_name}/stat.json'
        if input_type=='pe':
            print_input_type = 'Pair-End FASTQ'
        if input_type=='se':
            print_input_type = 'Single-End FASTQ'
        if input_type=='fa':
            print_input_type = 'FASTA'
        if os.path.isfile(stat_file) == False:
            worksheet.write(row, 0, sample_name, body_format)
            worksheet.write(row, 1, print_input_type, body_format)
            worksheet.write(row, 2, 'Skip', body_format)
            row += 1
            continue
        with open(stat_file) as handle:
            stat = json.load(handle)
            worksheet.write(row, 0, sample_name, body_format)
            worksheet.write(row, 1, print_input_type, body_format)
            # column = 2
            # for value in sorted(stat.values(), reverse=True):
            #     worksheet.write(row, column, value, body_format)
            #     column += 1
            worksheet.write(row, 2, stat['read_num_before_quality_filter'], body_format)
            worksheet.write(row, 3, stat['read_num_after_quality_filter'], body_format)
            worksheet.write(row, 4, stat['read_assemble'], body_format)
            worksheet.write(row, 5, stat['read_without_N'], body_format)
            worksheet.write(row, 6, stat['read_fit_length'], body_format)
            worksheet.write(row, 7, stat['uniq_read_num'], body_format)
            worksheet.write(row, 8, stat['haploid_num'], body_format)
            worksheet.write(row, 9, stat['species_num'], body_format)
        row += 1
    workbook.close()

def taxonomy(workdir:str, sample2type:dict):
    workbook = xlsxwriter.Workbook(f'{workdir}/MiFishResult/taxonomy.xlsx')
    with open(f'{workdir}/MiFishResult/species_meta.json') as in_handle:
        meta_data = json.load(in_handle)
    # Comparison between Samples
    header_format = workbook.add_format({'bg_color': 'yellow', 'font_name': 'Arial', 'align':'center', 'border': 1})
    body_format = workbook.add_format({'font_name': 'Arial', 'border': 1})
    worksheet = workbook.add_worksheet('Comparison of Samples')
    worksheet.set_column(0,2,12)
    worksheet.set_column(3,4,25)
    worksheet.set_column(5,5,15)
    worksheet.set_column(6,6,25)
    worksheet.set_column(7,7,15)
    worksheet.set_column(10,12,20)

    headers = 'Class	Order	Family	Scientific Name	Common Name	Ave. Confidence	Water area	Habitat	DepthS	DepthD	IUCN Red List Status	Importance in Fisheries	Threat to Humans'.split('\t')
    samples_list = list(sample2type.keys())
    headers += samples_list
    for col in range(len(headers)):
        worksheet.write(0, col, headers[col], header_format)
    classes_fish = 'Myxini Hyperoartia Chondrichthyes Cladistia Actinopteri'.split()
    classes_nonfish = 'Amphibia Lepidosauria Aves Mammalia unknown'.split()
    compare_data = {}
    species_name_data = {}
    for sample_name in samples_list:
        taxonomy_file = f'{workdir}/MiFishResult/Sample-{sample_name}/04_blast/{sample_name}.json'
        if os.path.isfile(taxonomy_file) == True:
            with open(taxonomy_file) as handle:
                taxonomy_data = json.load(handle)
                for (tax_name, tax_info) in taxonomy_data.items():
                    if tax_name not in compare_data:
                        compare_data[tax_name] = {}
                    compare_data[tax_name][sample_name] = {'total_read': tax_info['total_read'], 'main_confidence': tax_info['main_confidence']}
                    short_name = tax_info['short_name']
                    if short_name in meta_data:
                        classifications = meta_data[short_name]['classifications']
                        if len(classifications) > 0:
                            className = classifications[2]
                            order = classifications[3]
                            family = classifications[4]
                        else:
                            (className, order, family) = ('unknown', '', '')
                        if className not in species_name_data:
                            species_name_data[className] = {}
                        species_name_data[className][tax_name] = [tax_name, order, family, meta_data[short_name]['fishbase_data'], meta_data[short_name]['IUCN']]

    row = 1
    for className in classes_fish:
        if className in species_name_data:
            for tax_name in sorted(species_name_data[className].keys(), key = lambda j: (species_name_data[className][j][1], species_name_data[className][j][2])):
                (i, order, family, fishbase_data, IUCN_level) = species_name_data[className][tax_name]
                worksheet.write(row, 0, className, body_format)
                worksheet.write(row, 1, order, body_format)
                worksheet.write(row, 2, family, body_format)
                worksheet.write(row, 3, tax_name.replace('_', ' '), body_format)
                common_name = fishbase_data[1] if fishbase_data[1] != 'unknown' else ''
                worksheet.write(row, 4, common_name, body_format)

                abandance_data = compare_data[tax_name].values()
                average_confidence = max(abandance_data, key=lambda i:i['total_read'])['main_confidence']
                worksheet.write(row, 5, average_confidence, body_format)
                water = ''
                if fishbase_data[3] == 1:
                    water += 'Fresh Water; '
                if fishbase_data[5] == 1:
                    water += 'Salt Water; '
                if fishbase_data[4] == 1:
                    water += 'Brack Water; '
                worksheet.write(row, 6, water, body_format)
                habit = fishbase_data[6] if fishbase_data[6] != 'unknown' else ''
                worksheet.write(row, 7, habit, body_format)
                deep = fishbase_data[7] if fishbase_data[7] != 'unknown' else ''
                worksheet.write(row, 8, deep, body_format)
                deep = fishbase_data[8] if fishbase_data[8] != 'unknown' else ''
                worksheet.write(row, 9, deep, body_format)
                worksheet.write(row, 10, IUCN_level, body_format)
                importance = fishbase_data[9] if fishbase_data[9] != 'unknown' else ''
                worksheet.write(row, 11, importance, body_format)
                harm = fishbase_data[10] if fishbase_data[10] != 'unknown' else ''
                worksheet.write(row, 12, harm, body_format)
                col = 13
                for sample_name in samples_list:
                    size = compare_data[tax_name][sample_name]['total_read'] if sample_name in compare_data[tax_name] else 0
                    worksheet.write(row, col, size, body_format)
                    col += 1
                row += 1

    row += 2
    worksheet.write(row, 0, 'Followings are non-fish species', body_format)
    row += 1
    for className in classes_nonfish:
        if className in species_name_data:
            for tax_info in species_name_data[className]:
                (tax_name, order, family, fishbase_data, IUCN_raw_level) = species_name_data[className][tax_info]
                worksheet.write(row, 0, className, body_format)
                worksheet.write(row, 1, order, body_format)
                worksheet.write(row, 2, family, body_format)
                worksheet.write(row, 3, tax_name.replace('_', ' '), body_format)
                worksheet.write(row, 4, '', body_format)
                abandance_data = compare_data[tax_name].values()
                average_confidence = max(abandance_data, key=lambda i:i['total_read'])['main_confidence']
                worksheet.write(row, 5, average_confidence, body_format)
                worksheet.write(row, 6, '', body_format)
                worksheet.write(row, 7, '', body_format)
                worksheet.write(row, 8, '', body_format)
                worksheet.write(row, 9, '', body_format)
                worksheet.write(row, 10, '', body_format)
                worksheet.write(row, 11, '', body_format)
                worksheet.write(row, 12, '', body_format)
                col = 13
                for sample_name in samples_list:
                    size = compare_data[tax_name][sample_name]['total_read'] if sample_name in compare_data[tax_name] else 0
                    worksheet.write(row, col, size, body_format)
                    col += 1
                row += 1


    # List of Sample Details
    worksheet = workbook.add_worksheet('List of Sample Details')
    
    headers = 'Sample name	Species	Total read	Representative Sequence	Haploid ID	Size	Confidence	Identity(%)	LOD Score	Align Len	Mismatch	2nd-sp Name	2nd-sp Align Len	2nd-sp Mismatch	Sequence'.split('\t')
    worksheet.set_column(0,0,30)
    worksheet.set_column(1,1,25)
    worksheet.set_column(2,2,10)
    worksheet.set_column(3,3,30)
    worksheet.set_column(4,10,10)
    worksheet.set_column(11,11,25)
    worksheet.set_column(12,14,15)
    for col in range(len(headers)):
        worksheet.write(0, col, headers[col], header_format)
    
    body_merge_format = workbook.add_format({'font_name': 'Arial', 'align':'vcenter', 'border': 1})
    row = 0

    for sample_name in samples_list:
        for col in range(len(headers)):
            worksheet.write(row, col, headers[col], header_format)
        taxonomy_file = f'{workdir}/MiFishResult/Sample-{sample_name}/04_blast/{sample_name}.json'
        if os.path.isfile(taxonomy_file) == False:
            worksheet.write(row+1, 0, sample_name, body_format)
            worksheet.write(row+1, 1, 'Skip', body_format)
            row += 2
            continue
        with open(taxonomy_file) as handle:
            taxonomy_data = json.load(handle)
            for (tax_name, tax_info) in taxonomy_data.items():
                haploids = tax_info['hits']
                haploid_num = len(haploids)
                if haploid_num > 1:
                    worksheet.merge_range(row+1, 0, row+haploid_num, 0, sample_name, body_merge_format)
                    worksheet.merge_range(row+1, 1, row+haploid_num, 1, tax_name.replace('_', ' '), body_merge_format)
                    worksheet.merge_range(row+1, 2, row+haploid_num, 2, tax_info['total_read'], body_merge_format)
                    worksheet.merge_range(row+1, 3, row+haploid_num, 3, tax_info['representative_seq'], body_merge_format)
                else:
                    worksheet.write(row+1, 0, sample_name, body_format)
                    worksheet.write(row+1, 1, tax_name.replace('_', ' '), body_format)
                    worksheet.write(row+1, 2, tax_info['total_read'], body_format)
                    worksheet.write(row+1, 3, tax_info['representative_seq'], body_format)
                for haploid in haploids:
                    worksheet.write(row+1, 4, haploid['query'], body_format)
                    worksheet.write(row+1, 5, haploid['size'], body_format)
                    worksheet.write(row+1, 6, haploid['confidence'], body_format)
                    worksheet.write(row+1, 7, float('{0:.2f}'.format(haploid['top_hit']['identity']*100)), body_format)
                    if haploid['LOD_score'] != '':
                        worksheet.write(row+1, 8, float('{0:.2f}'.format(haploid['LOD_score'])), body_format)
                    else:
                        worksheet.write(row+1, 8, '', body_format)
                    worksheet.write(row+1, 9, haploid['top_hit']['aln_len'], body_format)
                    worksheet.write(row+1, 10, haploid['top_hit']['mismatch'], body_format)
                    if len(haploid['second_hit']) > 0:
                        worksheet.write(row+1, 11, haploid['second_hit']['species_name'].replace('_', ' '), body_format)
                        worksheet.write(row+1, 12, haploid['second_hit']['aln_len'], body_format)
                        worksheet.write(row+1, 13, haploid['second_hit']['mismatch'], body_format)
                    else:
                        worksheet.write(row+1, 11, '', body_format)
                        worksheet.write(row+1, 12, '', body_format)
                        worksheet.write(row+1, 13, '', body_format)
                    worksheet.write(row+1, 14, haploid['top_hit']['seq'], body_format)
                    row += 1
        row += 1
    
    # poor alignment
    worksheet = workbook.add_worksheet('Haploids with Low Identities')
    worksheet.set_column(0,0,30)
    worksheet.set_column(1,1,10)
    worksheet.set_column(2,2,25)
    worksheet.merge_range(0, 0, 0, 7, 'These Haploids are below the threshold of BLAST identity you defined', header_format)
    headers = 'Sample name	Haploid ID	Species	Total read	Identity(%)	Align Len	Mismatch	Sequence'.split('\t')
    for col in range(len(headers)):
        worksheet.write(2, col, headers[col], header_format)
    with open(f'{workdir}/MiFishResult/poor_alignment.json') as in_handle:
        row = 3
        for aln in json.load(in_handle):
            worksheet.write(row, 0, aln['sample'], body_format)
            worksheet.write(row, 1, aln['query'], body_format)
            worksheet.write(row, 2, aln['species_name'], body_format)
            worksheet.write(row, 3, aln['size'], body_format)
            worksheet.write(row, 4, aln['identity']*100, body_format)
            worksheet.write(row, 5, aln['aln_len'], body_format)
            worksheet.write(row, 6, aln['mismatch'], body_format)
            worksheet.write(row, 7, aln['seq'], body_format)
            row += 1
    workbook.close()

if __name__ == '__main__':
    sample2type = {'s1_':'pe', 's2_amplicon': 'se', 'error_':'pe', 's3_amplicon':'fa'}
    taxonomy(job_id='test', sample2type=sample2type)
