#!/usr/bin/env python3

#    host-snap - Removes specific viral sequences
#                from a host genome and generates a SNAP index
#
#    Copyright (C) 2017 Dillon Barker
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import subprocess
import sys
from functools import reduce
from io import BytesIO
from pathlib import Path
from typing import List, Iterator, NamedTuple

import pandas as pd
from Bio import SeqIO

Location = NamedTuple('Location',
                      (('contig', str), ('start', int), ('stop', int)))

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--cores', type=int, default=1,
                        help='Number of CPU cores to use [1]')

    parser.add_argument('--snap-seed', type=int, default=20,
                        help='Seed size for SNAP [20]')

    parser.add_argument('--identity', type=float, default=85,
                        help='Minimum percent identity for a match [85]')

    parser.add_argument('--evalue', type=float, default=1.0,
                        help='Maximum expect value for BLASTn [1.0]')

    parser.add_argument('--length', type=int, default=1000,
                        help='Minimum alignment length to subtract putative \
                              viral sequence from the host [1000]')

    parser.add_argument('viruses', nargs='+', type=Path,
                        help='Viruses (FASTA format) to subtract \
                              from the host sequence')

    parser.add_argument('host', type=Path, help='FASTA-formatted host genome')

    return parser.parse_args()

def makeblastdb(host: Path) -> None:
    '''Executes external program makeblastdb'''

    cmd = ('makeblastdb', '-in', str(host), '-dbtype', 'nucl')
    subprocess.run(cmd, check=True,
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)


def blast_search(virus: Path, host: Path,
                 evalue: float, cores: int) -> pd.DataFrame:

    headers = ('query', 'subject', 'identity', 'length', 'mismatch', 'gapopen',
               'query_start', 'query_end', 'subject_start', 'subject_end',
               'evalue', 'bitscore')

    blast = ('blastn', '-query', str(virus), '-db', str(host), '-outfmt', '10',
             '-evalue', str(evalue), '-num_threads', str(cores))

    output = subprocess.run(blast, stdout=subprocess.PIPE)

    df_output = pd.read_table(BytesIO(output.stdout), sep=',', header=None)

    df_output.columns = headers

    return df_output

def blast_viruses(viruses: List[Path], host: Path,
                  evalue: float, cores: int) -> pd.DataFrame:


    results = (blast_search(virus, host, evalue, cores) for virus in viruses)

    data = reduce(pd.DataFrame.append, results)

    return data

def parse_blast_output(length_threshold: int, max_evalue: float,
                       identity_threshold: float,
                       blast_output: pd.DataFrame) -> Iterator[Location]:

    for _, row in blast_output.iterrows():


        conditions = (row['identity'] >= identity_threshold,
                      row['evalue'] < max_evalue,
                      row['length'] >= length_threshold)

        if all(conditions):

            # 1-indexing -> 0-indexing
            loc = Location(contig=row['subject'],
                           start=row['subject_start'] - 1,
                           stop=row['subject_end'])

            yield loc

def host_mask(host_path: Path, length_threshold: int, max_evalue: float,
              identity_threshold: float, blast_output: pd.DataFrame) -> Path:

    def replace_sequence(sequence, start, stop):

        start_, stop_ = sorted((start, stop))

        span = stop_ - start_

        return sequence[:start_] + ('N' * span) + sequence[stop_:]

    with host_path.open('r') as host_fasta:

        host = SeqIO.to_dict(SeqIO.parse(host_fasta, 'fasta'))

    for location in parse_blast_output(length_threshold, max_evalue,
                                       identity_threshold, blast_output):

        seq = host[location.contig].seq

        host[location.contig].seq = replace_sequence(seq,
                                                     location.start,
                                                     location.stop)

    host_iterator = (host[contig] for contig in host)

    outname = host_path.with_suffix('.masked.fasta')

    with outname.open('w') as out:
        SeqIO.write(host_iterator, out, 'fasta')


    return outname

def snap(masked: Path, seed_size: int, cores: int) -> None:

    out = (masked.parent / masked.stem).with_suffix('.snap_index')
    cmd = ('snap-aligner', 'index', masked, out,
           '-s', seed_size, '-t{}'.format(cores))

    snap_cmd = [str(x) for x in cmd]

    try:
        subprocess.run(snap_cmd, check=True)

    except subprocess.CalledProcessError:

        print('Something blew up. Check SNAP output for hints.',
              file=sys.stderr)

        raise

def main():

    args = arguments()

    makeblastdb(args.host)

    blast_results = blast_viruses(args.viruses, args.host,
                                  args.evalue, args.cores)

    masked = host_mask(args.host, args.length, args.evalue,
                       args.identity, blast_results)

    snap(masked, args.snap_seed, args.cores)

if __name__ == '__main__':
    main()
