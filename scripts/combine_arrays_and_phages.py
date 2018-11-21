from collections import defaultdict

OUT_FILE = 'arrays_on_phages'

crispr_arrays = defaultdict(lambda: defaultdict(set))
for line in open("all_meta_info"):
    genome, array, contig, system = line.strip().split()
    crispr_arrays[genome][contig].add('_'.join([array, system]))

virsorter = defaultdict(lambda: defaultdict(set))
for line in open("phages.coords.tsv"):
    genome, contig, phage = line.strip().split(maxsplit=2)
    virsorter[genome][contig].add(phage)

out = open(OUT_FILE, 'w')
for genome, contig_arrays in crispr_arrays.items():
    contig_phages = virsorter[genome]
    for contig, arrays in contig_arrays.items():
        phages = contig_phages.get(contig)
        for array in arrays:
            _, _, coords, _, _ = array.split('_')
            s, e = map(int, coords.split('-', maxsplit=1))
            phage_id = "0"
            for phage in phages or []:
                if not phage:
                    continue
                elif len(phage) == 1:
                    phage_id = phage
                elif phage_id == "0":
                    phage_id, coords, _ = phage.split()
                    vs, ve = map(int, coords.split('-', maxsplit=1))
                    if not(vs <= s <= ve or vs <= e <= ve):
                        phage_id = "0"
                else:
                    continue
            print("%s\t%s\t%s" % (genome, array, phage_id), file=out)
out.close()
            
