from collections import defaultdict
import sys
from os import path
import argparse

#CRISPRMETA = "/home/hielke/bep/data/fp_crisprmeta/out_fp/"
#SPACERSFILE = "spacers.newnewfiltered.extra"
#EXTENSION = ".crisprdetect.fp"

parser = argparse.ArgumentParser(
        description="Program that removes spacers hitting arrays as defined in processed crispr meta output")
parser.add_argument("CRISPRMETAFOLDER", 
        help="Folder with processed crisprmeta output")
parser.add_argument("EXTENSION", 
        help="Optional extension for the crisprmeta output",
        default="", nargs="?",
        )
args = parser.parse_args()

def create_spacer_dict(genome):
    spacer_dict = defaultdict(set)
    file = path.join(args.CRISPRMETAFOLDER, genome + args.EXTENSION)
    if not path.isfile(file):
        return None
    for line in open(file):
        _, array, contig, _ = line.split()
        _, _, coords, direction = array.split('_')[-4:]
        s, e = map(int, coords.split('-'))
        if direction == 'R':
            s, e = e, s
        spacer_dict[contig].add((s, e))
    return spacer_dict
#sys.exit(1)
genome = ''
interesting = True
for line in sys.stdin:
    splitted_line = line.split()
    new_genome = splitted_line[0]
    if new_genome != genome:
        # we have a new genome
        genome = new_genome
        spacer_dict = create_spacer_dict(genome)
        interesting = spacer_dict is not None
    if not interesting:
        continue
    hc = int(splitted_line[2])
    contig = splitted_line[4]
    arrays = spacer_dict[contig]
    for arr in arrays:
        s, e = arr
        if s <= hc <= e:
            print(genome, file=sys.stderr)
            break
    else: # The loop executed normally
        print(line, end='')
