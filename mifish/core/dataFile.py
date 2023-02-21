from glob import glob
import os
import re
import sys
from typing import Tuple

def fetch(data_dirs:list, workdir:str)->Tuple[dict, dict]:
    name2file = {}
    group_to_sample = {}
    group_rank = 1
    sample2type = {}
    for data_dir in data_dirs:
        sample2name_this_group = {}
        group_name = f'Group{group_rank}'
        group_rank += 1
        for file in glob(f'{data_dir}/*'):
            if os.path.isdir(file) is True:
                continue
            if file=='':
                continue
            if file.endswith('.gz') is False and file.endswith('.bz2') is False and file.endswith('.xz') is False:
                continue
            name = os.path.basename(file)
            name2file[name] = file
            file_extension = os.path.splitext(name)[1]  # e.g: .gz, .bz2, .xz
            sample = re.sub(f'\D[12](_001)?\.f(ast)?q\{file_extension}$', '', name)
            sample = re.sub(f'\.f(ast)?q\{file_extension}$', '', sample)
            sample = re.sub(f'\.f(ast)?a\{file_extension}$', '', sample)
            sample = re.sub('\W+', '_', sample) # Remove all dangerous characters
            if sample not in sample2name_this_group:
                sample2name_this_group[sample] = []
            sample2name_this_group[sample].append(name)
        
        for (sample, names) in sample2name_this_group.items():
            if sample in sample2type:
                print(f'Error: duplicate sample name {sample} detected', file=sys.stderr)
                exit(1)
            names = list(set(names))
            if len(names)==2:
                (file1, file2) = (name2file[names[0]], name2file[names[1]])
                if (re.search(f'\D[1](_001)?\.f(ast)?q\.', file1) is not None \
                    and re.search(f'\D[2](_001)?\.f(ast)?q\.', file2) is not None):
                    file_extension = os.path.splitext(os.path.basename(file1))[1]
                    os.system(f'cp {file1} {workdir}/MiFishResult/data/{sample}.1.fq{file_extension}')
                    file_extension = os.path.splitext(os.path.basename(file2))[1]
                    os.system(f'cp {file2} {workdir}/MiFishResult/data/{sample}.2.fq{file_extension}')
                elif (re.search(f'\D[2](_001)?\.f(ast)?q\.', file1) is not None \
                        and re.search(f'\D[1](_001)?\.f(ast)?q\.', file2) is not None):
                    file_extension = os.path.splitext(os.path.basename(file2))[1]
                    os.system(f'cp {file2} {workdir}/MiFishResult/data/{sample}.1.fq{file_extension}')
                    file_extension = os.path.splitext(os.path.basename(file1))[1]
                    os.system(f'cp {file1} {workdir}/MiFishResult/data/{sample}.2.fq{file_extension}')
                sample2type[sample] = 'pe'
            else:
                file1 = str(name2file[names[0]])
                file_extension = os.path.splitext(os.path.basename(file1))[1]
                if file1.endswith(f'fasta{file_extension}') or file1.endswith(f'fa{file_extension}') or file1.endswith(f'fas{file_extension}'):
                    os.system(f'cp {file1} {workdir}/MiFishResult/data/{sample}.fa{file_extension}')
                    sample2type[sample] = 'fa'
                else:
                    os.system(f'cp {file1} {workdir}/MiFishResult/data/{sample}.fq{file_extension}')
                    sample2type[sample] = 'se'
            if group_name not in group_to_sample:
                group_to_sample[group_name] = []
            group_to_sample[group_name].append(sample)
    return (sample2type, group_to_sample)

if __name__ == '__main__':
    print(fetch(['dummy/data']))
