import json
import os

from collections import defaultdict

from skbio import diversity

def relative_abandance(workdir:str, group_to_sample:dict) -> dict:
    with open(f'{workdir}/MiFishResult/species_meta.json') as in_handle:
        meta_data = json.load(in_handle)
    
    result = {}
    unit_rank = 'Kingdom Phylum Class Order Family Genus'.split()
    all_samples = []
    for (group, samples) in group_to_sample.items():
        for sample_name in samples:
            workdir_sample = f'{workdir}/MiFishResult/Sample-{sample_name}'
            if os.path.isfile(f'{workdir_sample}/04_blast/{sample_name}.json') == False:
                continue
            with open(f'{workdir_sample}/04_blast/{sample_name}.json') as handle:
                abandance_data = json.load(handle)
                total_read_sample = sum([x['total_read'] for x in abandance_data.values()])
                for (tax_name, tax_info) in abandance_data.items():
                    total_read = tax_info['total_read']
                    species_name_short = '_'.join(tax_name.split('_')[0:2])
                    if species_name_short in meta_data:
                        taxonomy = meta_data[species_name_short]['classifications']
                        all_samples.append(f'{group}-{sample_name}')
                        for i in range(2, len(unit_rank)):
                            unit = unit_rank[i]
                            unit_name = taxonomy[i] if i < len(taxonomy) else 'Other'
                            if unit_name == '' or unit_name == ' ':
                                unit_name = 'Other'
                            if unit not in result:
                                result[unit] = {}
                            if unit_name not in result[unit]:
                                result[unit][unit_name] = {}
                            if f'{group}-{sample_name}' not in result[unit][unit_name]:
                                result[unit][unit_name][f'{group}-{sample_name}'] = 0
                            result[unit][unit_name][f'{group}-{sample_name}'] += total_read/total_read_sample
    all_samples = list(set(all_samples))
    for unit in result.keys():
        for unit_name in result[unit].keys():
            for sample_name in all_samples:
                if sample_name not in result[unit][unit_name]:
                    result[unit][unit_name][sample_name] = 0
    
    filtered_result = {}
    for (unit, unit_count) in result.items():
        rank = 0
        if unit not in filtered_result:
            filtered_result[unit] = {}
        for unit_name in sorted(unit_count.keys(), key=lambda i:max(unit_count[i].values()), reverse=True):
            rank += 1
            if rank <= 10:
                filtered_result[unit][unit_name] = result[unit][unit_name]
            else:
                if 'Other' not in filtered_result[unit]:
                    filtered_result[unit]['Other'] = {}
                for sample_name in result[unit][unit_name].keys():
                    if sample_name not in filtered_result[unit]['Other']:
                        filtered_result[unit]['Other'][sample_name] = 0
                    filtered_result[unit]['Other'][sample_name] += result[unit][unit_name][sample_name]
    
    return filtered_result


def eco_diversity(workdir, group_to_sample) -> dict:
    if len(group_to_sample) < 2:
        return {}
    all_samples = []
    all_species = []
    count_data = defaultdict(dict)
    diversity_result = defaultdict(dict)
    count_valid_group = 0
    for (group, samples) in group_to_sample.items():
        has_valid_samples_in_this_group = False
        for sample_name in samples:
            workdir_sample = f'{workdir}/MiFishResult/Sample-{sample_name}'
            if os.path.isfile(f'{workdir_sample}/04_blast/{sample_name}.json') is False:
                continue
            all_samples.append(sample_name)
            has_valid_samples_in_this_group = True
            with open(f'{workdir_sample}/04_blast/{sample_name}.json') as handle:
                abandance_data = json.load(handle)
                for (tax_name, tax_info) in abandance_data.items():
                    all_species.append(tax_name)
                    count_data[sample_name][tax_name] = tax_info['total_read']
        if has_valid_samples_in_this_group is True:
            count_valid_group += 1
    if count_valid_group < 2:
        return {}
    all_species = list(set(all_species))
    data_for_diversity = [[count_data[x][y] if y in count_data[x] else 0 for y in all_species] for x in all_samples]
    diversity_result['alpha']['simpson'] = dict(diversity.alpha_diversity('simpson', data_for_diversity, all_samples, validate=False))
    diversity_result['alpha']['chao1'] = dict(diversity.alpha_diversity('chao1', data_for_diversity, all_samples, validate=False))
    diversity_result['alpha']['shannon'] = dict(diversity.alpha_diversity('shannon', data_for_diversity, all_samples, validate=False))
    diversity_result['group'] = group_to_sample

    bc_dm = diversity.beta_diversity("braycurtis", data_for_diversity, all_samples)
    bc_dm_result = dict(bc_dm.to_series())
    bc_dm_result = {f'{key[0]}:{key[1]}':bc_dm_result[key] for key in bc_dm_result}
    diversity_result['beta'] = bc_dm_result
    return diversity_result