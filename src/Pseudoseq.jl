module Pseudoseq

using BioSequences, FASTX, Kmers, Random

#include("PuzzleMaker.jl")
#include("sequencing/Sequencing.jl")
include("Sequencing.jl")
export needed_sample_size, expected_coverage, X, x, bp, coverage_report, uncovered_positions, uncovered_regions
include("sequencing_views.jl")

include("Molecules.jl")
export Molecules, amplify,fragment, select, sequence, subsample, tag, flip
include("Reads.jl")
export paired_reads, unpaired_reads, nreads, ClearSubstitutions, FixedProbSubstitutions, edit_substitutions, make_substitutions, generate

include("DSL.jl")
export makereads

end # module
