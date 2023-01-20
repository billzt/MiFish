import os
import json
import sys
import math
import time
from glob import glob
import csv

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Blast import NCBIXML
import duckdb

import inspect

import mifish
from mifish.core import sizeFasIntegrator, excel, drawTree, stat, dataFile

genus2taxonomy = {}
with open(os.path.dirname(inspect.getfile(mifish)) + '/data/lineage.csv') as in_handle:
    csvReader = csv.DictReader(in_handle)
    for row in csvReader:
        genus2taxonomy[row['genus_name']] = [row['kingdom_name'], row['phylum_name'], row['class_name'], \
            row['order_name'], row['family_name'], row['genus_name']]

def calcConfidence(LOD_score, identity:float)->str:
    if LOD_score != '':
        if identity >= 0.98:
            if LOD_score >= 0.9:
                return 'HIGH'
            if 0.5 <= LOD_score < 0.9:
                return 'MODERATE'
            return 'LOW'
        else:
            return 'LOW'
    else:
        if identity >= 0.99:
            return 'HIGH'
        if 0.98 <= identity < 0.99:
            return 'MODERATE'
        return 'LOW'

def runMiFish(data_dir:str, data_dir_other_groups:list, min_read_length=204, max_read_length=254, \
    rm_p_5='GTCGGTAAAACTCGTGCCAGC', rm_p_3='CAAACTGGGATTAGATACCCCACTATG',\
        blast_identity=97.0, current_task='', current_message='', debug=True, simple_result=False, \
            workdir='.', cpu=2, db='mifishdb.fa'):

    all_species = {}
    low_hit_amplicons = []
    (sample2type, group_to_sample) = dataFile.fetch([data_dir] + data_dir_other_groups, workdir)
    print(f'Detect your data as', file=sys.stderr)
    for (group, samples) in group_to_sample.items():
        print(f'#########', file=sys.stderr)
        print(f'{group}: {len(samples)} samples', file=sys.stderr)
        for sample_name in samples:
            print(f'Sample {sample_name}: read type = {sample2type[sample_name]}', file=sys.stderr)
    print(f'#########', file=sys.stderr)
    for (sample_name, input_type) in sample2type.items():
        workdir_sample = f'{workdir}/MiFishResult/Sample-{sample_name}'
        stat_data = {}
        current_message = ''
        
        # process 0: files
        current_task = f'Sample {sample_name} Step 0: Decompress'
        if debug==True:
            print(current_task, file=sys.stderr)
        files = glob(f'{workdir}/MiFishResult/data/{sample_name}*')
        if input_type != 'fa':
            file_extension = os.path.splitext(os.path.basename(files[0]))[1]
            for file in files:
                if file_extension=='.bz2':
                    os.system(f'bzip2 -d {file}')
                if file_extension=='.xz':
                    os.system(f'xz -d {file}')
            files = glob(f'{workdir}/MiFishResult/data/{sample_name}*')

        # Step 1: filter the quality of FASTQ and merge Pair-End Reads
        current_task = f'Sample {sample_name} Step 1: filter the quality of FASTQ and merge Pair-End Reads'
        if debug==True:
            print(current_task, file=sys.stderr)
        if os.path.isdir(f'{workdir_sample}/01_filter_fastq_and_merge') is False:
            os.makedirs(f'{workdir_sample}/01_filter_fastq_and_merge')
        fastp_cmd = f'fastp --disable_quality_filtering --disable_length_filtering --disable_adapter_trimming\
            --dont_eval_duplication --thread 8 --report_title "Sample-{sample_name}" \
                --json {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.json \
                    --html {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.html \
                        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 10 \
                            --cut_tail --cut_front_window_size 1 --cut_tail_mean_quality 10'
        if input_type == 'pe':
            os.system(f'{fastp_cmd} --in1 {files[0]} --out1 \
                {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.1.fq.gz \
                    --in2 {files[1]} --out2 {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.2.fq.gz \
                        >{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.log 2>&1')
            os.system(f'flash {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.1.fq.gz {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.2.fq.gz \
                -d {workdir_sample}/01_filter_fastq_and_merge -o {sample_name} -z \
                    >{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.flash.log 2>&1')
            os.system(f'cat {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.extendedFrags.fastq.gz \
                {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.notCombined*.gz \
                >{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fq.gz')
        if input_type == 'se':
            os.system(f'{fastp_cmd} --in1 {files[0]} --out1 \
                {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fq.gz \
                    >{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.log 2>&1')
        if input_type != 'fa':
            os.system(f'seqkit fq2fa {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fq.gz \
                -o {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fa')
            # statics
            # make sure no error happens in fastp
            check_log = os.popen(f'file -b {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.log').read().strip()
            if check_log == 'data':
                current_message = f'Sample {sample_name} has error in FASTQ files. Skip'
                if debug==True:
                    print(current_message, file=sys.stderr)
                else:
                    time.sleep(1)
                continue
            with open(f'{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.log') as handle:
                lines = handle.readlines()
                if input_type == 'pe':
                    stat_data['read_num_before_quality_filter'] = int(lines[1].split(':')[1].strip())*2
                    stat_data['read_num_after_quality_filter'] = int(lines[13].split(':')[1].strip())*2
                else:
                    stat_data['read_num_before_quality_filter'] = int(lines[1].split(':')[1].strip())
                    stat_data['read_num_after_quality_filter'] = int(lines[7].split(':')[1].strip())
            # if input_type == 'pe':
            #     with open(f'{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.flash.log') as handle:
            #         assemble_percent = 0
            #         for line in handle:
            #             if 'Combined reads:' in line:
            #                 stat_data['read_assemble'] = int(line.split()[3].strip())*2
            #                 assemble_percent = stat_data['read_assemble']/stat_data['read_num_after_quality_filter']
            #                 break
            #         # if assemble_percent < 0.1:
            #         #     current_message = f"Sample {sample_name} has too few assembled \
            #         #         reads ({stat_data['read_assemble']}/{stat_data['read_num_after_quality_filter']}={assemble_percent:.1%}). Skip. \
            #         #             If the length of your pair-end read is close to or longer than the length of amplicons, please combine the pair-end file into \
            #         #                 a single one to treat them as single-end reads"
            #         #     if debug==True:
            #         #         print(current_message, file=sys.stderr)
            #         #     else:
            #         #         time.sleep(1)
            #         #     continue
            # else:
            #     stat_data['read_assemble'] = stat_data['read_num_after_quality_filter']
        else:
            file_extension = os.path.splitext(os.path.basename(files[0]))[1]
            os.system(f'cp {files[0]} \
                {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fa{file_extension}')
            if file_extension == '.gz':
                decompress_cmd = 'gzip'
            if file_extension == '.bz2':
                decompress_cmd = 'bzip2'
            if file_extension == '.xz':
                decompress_cmd = 'xz'
            os.system(f'{decompress_cmd} -d {workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fa{file_extension}')
            stat_data['read_num_before_quality_filter'] = 0
            stat_data['read_num_after_quality_filter'] = 0
        stat_data['read_assemble'] = 0

        
        # Step 2: filter read length and remove primers
        current_task = f'Sample {sample_name} Step 2: filter read length and remove primers'
        if debug==True:
            print(current_task, file=sys.stderr)
        stat_data['read_without_N'] = 0
        stat_data['read_fit_length'] = 0
        if os.path.isdir(f'{workdir_sample}/02_process_fasta') is False:
            os.makedirs(f'{workdir_sample}/02_process_fasta')
        with open(f'{workdir_sample}/01_filter_fastq_and_merge/{sample_name}.amplicon.fa') as handle:
            with open(f'{workdir_sample}/02_process_fasta/{sample_name}.filterlen.fa', mode='w') as out_handle:
                for (title, seq) in SimpleFastaParser(handle):
                    if input_type == 'fa':
                        stat_data['read_num_before_quality_filter'] += 1
                        stat_data['read_num_after_quality_filter'] += 1
                    stat_data['read_assemble'] += 1
                    if 'N' in seq:
                        continue
                    stat_data['read_without_N'] += 1
                    if min_read_length <= len(seq) <=  max_read_length:
                        stat_data['read_fit_length'] += 1
                        print(f'>{title}\n{seq}', file=out_handle)
        if stat_data['read_fit_length'] < 10:
            # This sample has not passed length filter
            current_message = f'Sample {sample_name} has not passed read length filter. Skip'
            if debug==True:
                print(current_message, file=sys.stderr)
            else:
                time.sleep(1)
            continue
        
        os.system(f'cutadapt -g "{rm_p_5};max_error_rate=0.15...{rm_p_3};max_error_rate=0.15" --revcomp -j {cpu} \
            {workdir_sample}/02_process_fasta/{sample_name}.filterlen.fa \
                >{workdir_sample}/02_process_fasta/{sample_name}.processed.fa \
                    2>{workdir_sample}/02_process_fasta/{sample_name}.cut.log')
        
        # Step 3: De-noise and generate haploid
        current_task = f'Sample {sample_name} Step 3: De-noise and generate haploid'
        if debug==True:
            print(current_task, file=sys.stderr)
        stat_data['uniq_read_num'] = 0
        stat_data['haploid_num'] = 0
        if os.path.isdir(f'{workdir_sample}/03_haploid') is False:
            os.makedirs(f'{workdir_sample}/03_haploid')
        os.system(f'usearch -fastx_uniques {workdir_sample}/02_process_fasta/{sample_name}.processed.fa \
            -sizeout -relabel Uniq -fastaout {workdir_sample}/03_haploid/{sample_name}.derep.fasta -threads {cpu} \
                >{workdir_sample}/03_haploid/{sample_name}.fastx_uniques.log 2>&1')
        os.system(f'usearch -unoise3 {workdir_sample}/03_haploid/{sample_name}.derep.fasta \
            -zotus {workdir_sample}/03_haploid/{sample_name}.zotus.fasta -threads {cpu} \
                -tabbedout {workdir_sample}/03_haploid/{sample_name}.zotus.details.txt \
                    >{workdir_sample}/03_haploid/{sample_name}.unoise3.log 2>&1')
        if os.path.isfile(f'{workdir_sample}/03_haploid/{sample_name}.zotus.details.txt') is None \
            or os.path.getsize(f'{workdir_sample}/03_haploid/{sample_name}.zotus.details.txt')==0:
            # This sample has no haploids
            current_message = f'Sample {sample_name} has no haploids. Skip'
            if debug==True:
                print(current_message, file=sys.stderr)
            else:
                time.sleep(1)
            continue
        os.system(f'usearch -otutab {workdir_sample}/02_process_fasta/{sample_name}.processed.fa \
            -zotus {workdir_sample}/03_haploid/{sample_name}.zotus.fasta -threads {cpu} \
                -otutabout {workdir_sample}/03_haploid/{sample_name}.zotus.size.txt \
                    >{workdir_sample}/03_haploid/{sample_name}.otutab.log 2>&1')
        sizeFasIntegrator.run(zotusCountFile=f'{workdir_sample}/03_haploid/{sample_name}.zotus.size.txt', \
            zotusFaFile=f'{workdir_sample}/03_haploid/{sample_name}.zotus.fasta', \
                resultFile=f'{workdir_sample}/03_haploid/{sample_name}.zotus.size.fasta')
        os.system(f'usearch -sortbysize {workdir_sample}/03_haploid/{sample_name}.zotus.size.fasta \
            -fastaout {workdir_sample}/03_haploid/{sample_name}.sortsize.fasta -threads {cpu} \
                >{workdir_sample}/03_haploid/{sample_name}.sortbysize.log 2>&1')
        stat_data['uniq_read_num'] = int(os.popen(f'grep -c ">" {workdir_sample}/03_haploid/{sample_name}.derep.fasta').read().strip())
        stat_data['haploid_num'] = int(os.popen(f'grep -c ">" {workdir_sample}/03_haploid/{sample_name}.sortsize.fasta').read().strip())

        # Step 4: BLAST and calculate LOD Score
        current_task = f'Sample {sample_name} Step 4: BLAST and calculate LOD Score'
        if debug==True:
            print(current_task, file=sys.stderr)
        if os.path.isdir(f'{workdir_sample}/04_blast') is False:
            os.makedirs(f'{workdir_sample}/04_blast')
        os.system(f'blastn -num_threads {cpu} -query {workdir_sample}/03_haploid/{sample_name}.sortsize.fasta -db {db} \
            -max_target_seqs 5 -evalue 0.00001 -outfmt 5 \
                -out {workdir_sample}/04_blast/{sample_name}.xml')
        query_seq_obj = {}
        with open(f'{workdir_sample}/03_haploid/{sample_name}.sortsize.fasta') as handle:
            for (title, seq) in SimpleFastaParser(handle):
                query_seq_obj[title] = seq
        species2amplicon = {}
        with open(f'{workdir_sample}/04_blast/{sample_name}.xml') as handle:
            for blast_record in NCBIXML.parse(handle):
                query_seq = query_seq_obj[blast_record.query+';']
                (query, size) = blast_record.query.split(';')
                size = int(size.replace('size=', ''))
                top_hit = {}
                second_hit = {}
                LOD_score = ''
                good_alns = []
                poor_alns = []
                for alignment in blast_record.alignments:
                    hsp = alignment.hsps[0]
                    aln_len = alignment.length
                    identity = hsp.identities/aln_len
                    if identity >= blast_identity/100:
                        good_alns.append(alignment)
                    else:
                        poor_alns.append(alignment)
                        break
                for alignment in good_alns:
                    aln_len = alignment.length
                    hit = alignment.hit_def
                    species_name = hit.split('|')[2].replace(':', '_')
                    hsp = alignment.hsps[0]
                    identity = hsp.identities/aln_len
                    mismatch = aln_len - hsp.identities
                    if len(top_hit) == 0:
                        top_hit = {'species_name':species_name, 'aln_len':aln_len, \
                            'identity':identity, 'mismatch':mismatch, 'seq':query_seq}
                    elif top_hit['species_name'] != species_name:
                        second_hit = {'species_name':species_name, 'aln_len':aln_len, \
                            'identity':identity, 'mismatch':mismatch, 'seq':query_seq}
                        break
                if len(second_hit) > 0:
                    LOD_score = math.log((top_hit['aln_len']/(top_hit['mismatch']+1))/(second_hit['aln_len']/(second_hit['mismatch']+1)))
                if len(top_hit) > 0:
                    if top_hit['species_name'] not in species2amplicon:
                        species2amplicon[top_hit['species_name']] = []
                    species2amplicon[top_hit['species_name']].append({'top_hit':top_hit, \
                        'second_hit':second_hit, 'LOD_score':LOD_score, \
                            'confidence':calcConfidence(LOD_score, top_hit['identity']), 'size':size})
                if len(good_alns) == 0:
                    for alignment in poor_alns:
                        aln_len = alignment.length
                        hit = alignment.hit_def
                        species_name = hit.split('|')[2].replace(':', '_')
                        hsp = alignment.hsps[0]
                        identity = hsp.identities/aln_len
                        mismatch = aln_len - hsp.identities
                        low_hit_amplicons.append({'sample': sample_name, 'query':query, 'species_name':species_name, 'aln_len':aln_len, \
                                'identity':identity, 'mismatch':mismatch, 'size': size, 'seq':query_seq})
        species_result = {}
        for (species_name, top_hits) in species2amplicon.items():
            total_read = 0
            confidence_to_size = {'HIGH':0, 'MODERATE':0, 'LOW':0}
            representative_seq = ''
            representative_size = 0
            for x in top_hits:
                total_read += x['size']
                confidence_to_size[x['confidence']] += x['size']
                if x['size'] > representative_size:
                    representative_size = x['size']
                    representative_seq = x['top_hit']['seq']
            main_confidence = max(confidence_to_size.keys(), key=lambda x:confidence_to_size[x])
            species_result[species_name] = {'short_name':'_'.join(species_name.split('_')[0:2]), 'total_read':total_read, 'taxonomy':'', 'JPname': '', \
                'main_confidence':main_confidence, 'hits':top_hits, 'representative_seq':representative_seq}
            species_name_short = '_'.join(species_name.split('_')[0:2])
            all_species[species_name_short] = 1
        with open(f'{workdir_sample}/04_blast/{sample_name}.json', mode='w') as out_handle:
            json.dump(species_result, fp=out_handle, indent=4, ensure_ascii=False)
        stat_data['species_num'] = len(species_result)
        if stat_data['species_num'] == 0:
            # This sample has no fish
            current_message = f'Sample {sample_name} has no fishes or other species in DB. Skip'
            if debug==True:
                print(current_message, file=sys.stderr)
            else:
                time.sleep(1)
            continue

        # Step 5: Phylogenetic Analysis
        if simple_result == False and stat_data['species_num'] > 3:
            current_task = f'Sample {sample_name} Step 5: Phylogenetic Analysis'
            if debug==True:
                print(current_task, file=sys.stderr)
            if os.path.isdir(f'{workdir_sample}/05_MSA') is False:
                os.makedirs(f'{workdir_sample}/05_MSA')
            with open(f'{workdir_sample}/05_MSA/{sample_name}.fa', mode='w') as out_handle:
                for (species_name, species_info) in species_result.items():
                    amplicon_rank = 1
                    hits = species_info['hits']
                    for hit in hits:
                        if len(hits)==1:
                            print(f'>{hit["top_hit"]["species_name"]}\n{hit["top_hit"]["seq"]}', file=out_handle)
                        else:
                            print(f'>{hit["top_hit"]["species_name"]}[{amplicon_rank}]\n{hit["top_hit"]["seq"]}', file=out_handle)
                        amplicon_rank += 1
            os.system(f'mafft --thread -1 --quiet --preservecase --adjustdirection \
                {workdir_sample}/05_MSA/{sample_name}.fa >{workdir_sample}/05_MSA/{sample_name}.aln.fa 2>{workdir_sample}/05_MSA/mafft.log')
            os.system(f'Gblocks {workdir_sample}/05_MSA/{sample_name}.aln.fa -t=d -b4=5 -b5=h -e=.gb \
                >{workdir_sample}/05_MSA/gblock.log 2>&1')
            os.system(f'FastTreeMP -gtr -nt {workdir_sample}/05_MSA/{sample_name}.aln.fa.gb \
                >{workdir_sample}/05_MSA/{sample_name}.nwk 2>{workdir_sample}/05_MSA/FastTree.log')
            drawTree.svg(species_result=species_result, tree_file=f'{workdir_sample}/05_MSA/{sample_name}.nwk', \
                svg_file=f'{workdir_sample}/05_MSA/{sample_name}.svg', sample_name=sample_name)

        with open(f'{workdir_sample}/stat.json', mode='w') as out_handle:
            json.dump(stat_data, fp=out_handle, indent=4)
        
    
    # Species Meta Data
    conn = duckdb.connect()
    fishbase_file = os.path.dirname(inspect.getfile(mifish)) + '/data/species.parquet'
    stock_file = os.path.dirname(inspect.getfile(mifish)) + '/data/stocks.parquet'
    link_fishbase = conn.from_parquet(fishbase_file)
    link_stock = conn.from_parquet(stock_file)
    IUCN_levels = {'N.E.': 'Not Evaluated', 'DD': 'Data Deficient', 'N.A.': 'Not Available', \
        'LC': 'Least Concern', 'LR/lc': 'Lower Risk: least concern', 'NT': 'Near Threatened', \
            'LR/cd': 'Vulnerable', 'VU': 'Vulnerable', 'LR/nt': 'Lower Risk: near threatened', \
                'EN': 'Endangered', 'EX': 'Extinct', 'EW': 'Extinct in the Wild', 'CR': 'Critically Endangered'}
    for species_name in all_species.keys():
        (genus_name, *species_subnames) = species_name.split('_')
        species_subname = species_subnames[0] if len(species_subnames) > 0 else 'sp'
        fishbase_data = link_fishbase.filter(f"Genus = '{genus_name}' AND Species = '{species_subname}'").project('SpecCode, FBname, \
                PicPreferredName, Fresh, Brack, Saltwater, DemersPelag, \
                    DepthRangeShallow, DepthRangeDeep, Importance, Dangerous').fetchone()
        if fishbase_data is not None:
            fishbase_data = list(fishbase_data)
            for i in range(len(fishbase_data)):
                if fishbase_data[i] is None:
                    fishbase_data[i] = 'unknown'
            IUCN_data = link_stock.filter(f"SpecCode = {fishbase_data[0]}").project('IUCN_Code').fetchone()
            IUCN_data_level = IUCN_levels[IUCN_data[0]] if IUCN_data[0] in IUCN_levels else 'Not Available'
        else:
            fishbase_data = ['unknown']*11
            IUCN_data_level = 'Not Available'
        classifications = genus2taxonomy[genus_name] if genus_name in genus2taxonomy else []
        all_species[species_name] = {'fishbase_data': fishbase_data, 'IUCN': IUCN_data_level, 'classifications':classifications}

    conn.close()
    with open(f'{workdir}/MiFishResult/species_meta.json', mode='w') as out_handle:
        json.dump(all_species, fp=out_handle, indent=4)

    # Relative Abandance
    if simple_result == False:
        with open(f'{workdir}/MiFishResult/relative_abandance.json', mode='w') as out_handle:
            json.dump(stat.relative_abandance(workdir, group_to_sample), fp=out_handle, indent=4)
    
    # Diversity
    if simple_result == False:
        with open(f'{workdir}/MiFishResult/diversity.json', mode='w') as out_handle:
            json.dump(stat.eco_diversity(workdir, group_to_sample), fp=out_handle, indent=4)

    # poor amplicons
    with open(f'{workdir}/MiFishResult/poor_alignment.json', mode='w') as out_handle:
        json.dump(low_hit_amplicons, fp=out_handle, indent=4)

    # Excel
    excel.read_stat(workdir, sample2type)
    excel.taxonomy(workdir, sample2type)

    # QC HTML, NWK Tree
    os.system(f'zip -j {workdir}/MiFishResult/QC.zip {workdir}/MiFishResult/Sample-*/01_filter_fastq_and_merge/*.html')
    if simple_result == False and stat_data['species_num'] > 3:
        os.system(f'zip -j {workdir}/MiFishResult/tree.zip {workdir}/MiFishResult/Sample-*/05_MSA/*.nwk {workdir}/MiFishResult/Sample-*/05_MSA/*.svg')

    current_task = 'Finished'

    # Remove unnecessary files
    os.system(f'rm -rf {workdir}/MiFishResult/poor_alignment.json {workdir}/MiFishResult/species_meta.json {workdir}/MiFishResult/Sample-*')