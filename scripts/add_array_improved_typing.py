from collections import defaultdict 
from operator import itemgetter

new_arrays = defaultdict(lambda: defaultdict(list))
for l in open("out_crispr_new_meta"):
    genome, array, contig, repeat_type = l.strip('\n').split('\t')
    new_arrays[genome][contig].append("_".join([array, repeat_type]))

arrays_to_types = defaultdict(lambda: defaultdict(list))
with open("re_systems_comb.6.m") as f:
    next(f) # header 
    for l in f:
        F = l.rstrip('\n').split('\t')
        genome = F[0]
        conf = F[-1]
        arrays = F[-2]
        crispr_type = F[10]
        if not crispr_type:
            crispr_type = F[1]
        additional_information = ",".join(F[2:5])
        for a, c in zip(arrays.split(','), conf.split(',')):
            if a:
                if c == "3": # we might be able to improve the classification
                    crispr_types = crispr_type.split(",")
                    if len(crispr_types) != 1: # now we can
                        repeat_system = a.split("_")[-1]
                        repeat_systems = ['-'.join([repeat_system.split('-')[0], p]) for p in repeat_system.split('-')[1].split('/')]
                        res = [t for t in crispr_types if t[4:] in repeat_systems]
                        crispr_type = ",".join(res)
                #if a == "Array_1_1571463-1571190_R_II-C":
                #    raise Exception
                arrays_to_types[genome][a].append((c, crispr_type, additional_information))

prev_genome = ''
for l in open("all_meta_info_enhanced"):
    l = l.rstrip('\n')
    F = l.split("\t")
    genome, _, contig = F[:3]
    arrays_to_types_on_genomes = arrays_to_types_on_genomes if \
            prev_genome == genome else arrays_to_types[genome]
    c_1, c_2 = map(int, F[-2:])
    new_arrays_on_genomes = new_arrays_on_genomes if \
            prev_genome == genome else new_arrays[genome]
    prev_genome = genome
    for a in new_arrays_on_genomes[contig]:
        a_1, a_2 = map(int, a.split("_")[2].split("-"))
        # If the total width of the range containing the two ranges is smaller
        # than the sum of the two ranges, then the ranges overlap.
        overlaid = max([a_1, a_2, c_1, c_2]) - min([a_1, a_2, c_1, c_2]) < \
                abs(a_1 - a_2) + abs(c_1 - c_2)
        if not overlaid:
            continue
        lst = arrays_to_types_on_genomes[a]
        if lst:
            if len(lst) == 1:
                conf = lst[0][0]
                crispr_types = [lst[0][1]]
                additional_informations = [lst[0][2]]
            else:
                conf = max(lst, key=itemgetter(0))[0]
                crispr_types = [i[1] for i in lst if i[0] == conf]
                additional_informations = [i[2] for i in lst if i[0] == conf]
        else:
            conf = ""
            crispr_types = [""]
            additional_informations = [""]
        print("\t".join([l, a.split("_")[-1], a, conf, "&".join(additional_informations), "&".join(crispr_types)]))
        break
