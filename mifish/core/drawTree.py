from ete3 import Tree, TreeStyle, TextFace, NodeStyle, CircleFace
import os
import re
os.environ['QT_QPA_PLATFORM']='offscreen'

def svg(species_result:dict, tree_file:str, svg_file:str, sample_name:str):
    total_size = sum([x['total_read'] for x in species_result.values()])
    species_relative_abandance = {}
    unit_to_species = {}
    for (species_name, species_info) in species_result.items():
        amplicon_rank = 1
        unit = species_info['taxonomy'][3] if len(species_info['taxonomy'])>0 else 'Other'
        if unit == '':
            unit = 'Other'
        if unit not in unit_to_species:
            unit_to_species[unit] = []
        hits = species_info['hits']
        for hit in hits:
            size = hit['size']
            percent = size/total_size
            if len(hits)==1:
                species_relative_abandance[species_name] = percent
                unit_to_species[unit].append(species_name)
            else:
                species_relative_abandance[f'{species_name}[{amplicon_rank}]'] = percent
                unit_to_species[unit].append(f'{species_name}[{amplicon_rank}]')
            amplicon_rank += 1
    
    tree_handle = Tree(tree_file)
    for node in tree_handle.traverse():
        if node.is_leaf() == False:
            if node.support < 0.7:
                sp_value = ''
            else:
                sp_value = '{0:.0f}'.format(node.support*100)
            node.add_face(TextFace(sp_value, ftype='Arial', fsize=3), column=0, position='branch-top')
            node.set_style(NodeStyle(size=1, fgcolor='black', hz_line_width=1, vt_line_width=1))
        else:
            clean_name = re.sub('^_R_', '', node.name)
            node.name = clean_name
            i = TextFace(node.name.replace('_', ' '), ftype='Arial', fsize=5)
            i.margin_left = 5
            node.add_face(i, column=0)
            percent = species_relative_abandance[node.name]
            if percent>0.1:
                percent_color = 'red'
                percent_size = 8
            elif percent>0.05:
                percent_color = 'orange'
                percent_size = 6
            elif percent>0.01:
                percent_color = 'green'
                percent_size = 4
            else:
                percent_color = 'gray'
                percent_size = 3
            node.set_style(NodeStyle(fgcolor=percent_color, size=percent_size, hz_line_width=1, vt_line_width=1))

    for (unit, leaf_names) in unit_to_species.items():
        node = tree_handle.get_common_ancestor(leaf_names)
        if node.is_root() == False:
            i = TextFace(unit, ftype='Arial', fgcolor='red', fsize=3)
            i.margin_top = 5
            i.opacity = 0.5
            node.add_face(i, column=0, position='float-behind')

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.legend_position = 2
    ts.margin_left = 5
    ts.margin_right = 15
    ts.margin_top = 5
    ts.margin_bottom = 5
    ts.title.add_face(TextFace(f"Sample: {sample_name}", ftype='Arial', fsize=4), column=0)
    ts.legend.add_face(TextFace("Relative ", ftype='Arial', fsize=4), column=0)
    ts.legend.add_face(TextFace("Abandance", ftype='Arial', fsize=4), column=1)
    ts.legend.add_face(CircleFace(radius=4, color='red'), column=0)
    ts.legend.add_face(TextFace(">10%", ftype='Arial', fsize=4), column=1)
    ts.legend.add_face(CircleFace(radius=3, color='orange'), column=0)
    ts.legend.add_face(TextFace("5%~10%", ftype='Arial', fsize=4), column=1)
    ts.legend.add_face(CircleFace(radius=2, color='green'), column=0)
    ts.legend.add_face(TextFace("1%~5%", ftype='Arial', fsize=4), column=1)
    ts.legend.add_face(CircleFace(radius=1.5, color='gray'), column=0)
    ts.legend.add_face(TextFace("<1%", ftype='Arial', fsize=4), column=1)
    ts.show_scale = False
    tree_handle.render(svg_file, tree_style=ts, w=600, units="px")