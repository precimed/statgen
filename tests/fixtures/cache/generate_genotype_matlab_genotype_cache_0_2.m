% Generate tests/fixtures/cache/genotype_matlab_genotype_cache_0_2.mat.
% Uses the shared allele-swap dataset with indel swaps enabled.
% Fixture writer: GNU Octave 8.4.0, statgen MAT-file format v7.
% Run from the repository root with Octave (MAT-file format v7):
%   octave --no-gui --quiet tests/fixtures/cache/generate_genotype_matlab_genotype_cache_0_2.m

addpath(fileparts(mfilename('fullpath')));
cache_fixture_helpers('genotype');
