# » awk '$NF>75{match($1,/\/(.*).fna/,m); print m[1]"\t"$2"\t"$(NF-6)"\t"$(NF-5)"\t"$3"\t"$(NF-2)"\t"$NF}' all_spacers > spacers.filtered.with.hit

from collections import defaultdict

no_virsorter = None
with open('no_virsorter') as f:
    no_virsorter = set(map(lambda s: s.strip(), f))

with_phages = None
with open('genomes_with_phages') as f:
    with_phages = set(map(lambda s: s.strip(), f))

spacercount = defaultdict(dict)
with open('spacercount') as f:
    for line in f:
        genome, spacer_array, spacer_count = line.split()
        spacercount[genome][spacer_array] = int(spacer_count)

phages_dict = defaultdict(lambda : defaultdict(list))
with open('phages.coords.tsv') as f:
    for line in f:
        genome, contig, what = line.split(maxsplit=2)
        phages_dict[genome][contig].append(what.strip())

new_file = open('spacers.filtered.extra', 'w')
with open('spacers.filtered.with.hit') as f:
    phages_what = None
    spacercount_genome = None
    old_genome = None
    for line in f:
        genome, _, c_1, c_2, hit_contig, _, _ = line.split()
        if old_genome != genome:
            phages_what = phages_dict[genome]
            spacercount_genome = spacercount[genome]
            old_genome = genome
        if not phages_what:
            if genome in no_virsorter:
                phage_id = '-1'
            else:
                phage_id = '0'
        else:
            phages = phages_what[hit_contig]
            if not phages:
                phage_id = "0"
            else:
                for phage in phages:
                    if len(phage) == 1:
                        phage_id = phage
                    else:
                        phage_id, coords, _ = phage.split()
                        start, end = map(int, coords.split('-'))
                        c_1, c_2 = int(c_1), int(c_2)
                        if not(start <= c_1 <= end or start <= c_2 <= end):
                            phage_id = '0'
                        else: # its got its signature
                            break

        with_phage = '1' if genome in with_phages else phage_id 
        #phage_id is '0' or '-1' if there are no phages: 
        #it is not there, or unknown, resp.

        this_spacer_count = sum(spacercount_genome.values())

        new_line = '%s\t%d\t%s\t%s\n' % (line.strip(), this_spacer_count,
                phage_id, with_phage)
        new_file.write(new_line)

new_file.close()
