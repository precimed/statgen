"""Run from repository root AT source commit c36657c35d762bd6a92f1175313bac4a51167088.

Captures historical kernels/caches before implementation. New expected hashes
come only from the scalar recurrence below, never from production kernels.
Existing historical artifacts must not be regenerated with newer writers.
"""
import csv
from pathlib import Path
import subprocess
import warnings

from statgen._utils import allele_hash64, _ALLELE_HASH_MAX_CHARS
from statgen.reference import load_reference, save_reference_cache

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def scalar_hash(allele):
    h1 = h2 = 1
    for byte in allele[:1000000].encode('ascii'):
        h1 = (h1 * 257 + byte + 1) % 2147483647
        h2 = (h2 * 263 + byte + 1) % 2147483647
    return (h1 << 32) + h2


def sequence(length, motif):
    return (motif * ((length + len(motif) - 1) // len(motif)))[:length]


def main():
    if _ALLELE_HASH_MAX_CHARS != 150:
        raise RuntimeError("Historical artifacts require the recorded source commit, with the 150-base kernels")
    cases = [(1, base) for base in 'ACGT']
    cases += [(n, 'ACGTGCA') for n in list(range(2, 152)) + [8191, 8192, 8193, 20000, 100000, 193000, 999999, 1000000, 1000001]]
    alleles = [sequence(n, motif) for n, motif in cases]
    (OUT / 'alleles.txt').write_text('\n'.join(alleles) + '\n')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        old = allele_hash64(alleles)
    with (OUT / 'vectors.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['length', 'motif', 'historical', 'expected'])
        writer.writerows((n, motif, int(h), scalar_hash(s)) for (n, motif), h, s in zip(cases, old, alleles))
    bim = OUT / 'legacy.bim'
    bim.write_text('1\tlegacy\t0\t100\t' + sequence(193000, 'ACGTGCA') + '\tT\n')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        save_reference_cache(load_reference(bim), OUT / 'legacy_python.npz')
    r = '''for (f in list.files("R-package/R", pattern="[.]R$", full.names=TRUE)) source(f)
    d <- "tests/fixtures/long_alleles/"
    a <- readLines(paste0(d,"alleles.txt"))
    writeLines(as.character(suppressWarnings(.allele_hash64(a))), paste0(d,"historical_r.txt"))
    save_reference_cache(suppressWarnings(load_reference(paste0(d,"legacy.bim"))),paste0(d,"legacy_r.rds"))
    '''
    subprocess.run(['Rscript', '-e', r], cwd=ROOT, check=True)
    m = '''addpath('matlab'); d='tests/fixtures/long_alleles/';
    fid=fopen([d 'alleles.txt']); a=textscan(fid,'%s'); fclose(fid);
    h=statgen.internal.allele_hash64(a{1}); fid=fopen([d 'historical_octave.txt'],'w'); fprintf(fid,'%u\\n',h); fclose(fid);
    r=statgen.load_reference([d 'legacy.bim']); statgen.save_reference_cache(r,[d 'legacy_matlab.mat'],'format','v7');
    '''
    subprocess.run(['octave', '--no-gui', '--quiet', '--eval', m], cwd=ROOT, check=True)
    for runtime in ['r', 'octave']:
        path = OUT / f'historical_{runtime}.txt'
        assert [int(x) for x in path.read_text().splitlines()] == [int(x) for x in old]
        path.unlink()
    (OUT / 'alleles.txt').unlink()


if __name__ == '__main__':
    main()
