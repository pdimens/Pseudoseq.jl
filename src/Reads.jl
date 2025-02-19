# Paired end reads
# ================

abstract type Pairing end
struct Paired <: Pairing end
struct UnPaired <: Pairing end

struct Substitution
    pos::UInt64
    base::DNA
end

struct Reads{P<:Pairing,V<:AbstractSequencingView}
    genome::Vector{LongSequence{DNAAlphabet{2}}}
    views::Vector{V}
    substitutions::Dict{Int,Vector{Substitution}}
end

@inline views(x::Reads) = x.views
@inline genome(x::Reads) = x.genome
@inline substitutions(x::Reads) = x.substitutions
@inline nreads(x::Reads) = length(x.views)
@inline pick_read(reads::Reads) = rand(1:nreads(reads))
@inline function nsubstitutions(x::Reads)
    subs = substitutions(x)
    return isempty(subs) ? 0 : sum(length(x) for x in values(subs))
end

# I/O and printing functions
# --------------------------

summary(x::Reads{Paired,T}) where {T} = string(nreads(x), " paired reads:")
summary(x::Reads{UnPaired,T}) where {T} = string(nreads(x), " unpaired reads:")

function printlen(x::Reads)
    mx, av, mn = summarize_lengths(views(x))
    return string(" Maximum read size: $mx\n Average read size: $av\n Minimum read size: $mn")
end

function Base.show(io::IO, reads::Reads)
    println(summary(reads))
    println(printlen(reads))
    println(string(" Number of errors: ", nsubstitutions(reads)))
end

# Constructors
# ------------

function Reads{Paired,V}(p::Molecules{V}, flen::Int, rlen::Int = flen) where {V<:AbstractSequencingView}
    vs = views(p)
    vse = take_paired_ends(vs, flen, rlen)
    rds = Reads{Paired,V}(genome(p), vse, Dict{Int,Vector{Substitution}}()) 
    return rds
end

function Reads{UnPaired,V}(p::Molecules{V}, len::Union{Int,Nothing}) where {V<:AbstractSequencingView}
    vs = views(p)
    vse = take_single_ends(vs, len)
    rds = Reads{UnPaired,V}(genome(p), vse, Dict{Int, Vector{Substitution}}())
    return rds
end

"""
    paired_reads(p::MoleculePool, flen::Int, rlen::Int = flen)

Create a set of paired-end reads from a pool of DNA molecules `p`.

`flen` sets the length of forward read, and `rlen` sets the length of the
reverse read. If you only provide `flen`, then the function sets `rlen = flen`.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function paired_reads(p::Molecules{T}, flen::Int, rlen::Int = flen) where {T<:AbstractSequencingView}
    return Reads{Paired,T}(p, flen, rlen)
end

"""
    unpaired_reads(p::Molecules{T}, len::Int) where {T<:AbstractSequencingView}

Create a set of single-end reads from a pool of DNA molecules `p`.

`len` sets the length of the reads.

The end (strand) from which the reading begins for each DNA molecule in the
pool is determined at random for each molecule, with 50:50 probability.

If you don't provide a value for `len`, then the function will read each
DNA molecule in it's entirety.

!!! note
    If a molecule in the pool is not long enough to create a forward and/or
    reverse read, then that molecule will simply be skipped. 
"""
function unpaired_reads(p::Molecules{T}, len::Union{Int,Nothing}) where {T<:AbstractSequencingView}
    return Reads{UnPaired,T}(p, len)
end

const posbases = (
    (DNA_C, DNA_G, DNA_T),
    (DNA_A, DNA_G, DNA_T),
    (DNA_N, DNA_N, DNA_N),
    (DNA_A, DNA_C, DNA_T),
    (DNA_N, DNA_N, DNA_N),
    (DNA_N, DNA_N, DNA_N),
    (DNA_N, DNA_N, DNA_N),
    (DNA_A, DNA_C, DNA_G)
)

"""
A simple function object that can be used with [`edit_substitutions`](@ref) or
[`edit_substitutions!`](@ref).

Simply clears all sequencing errors - substitutions - for all reads.
"""
struct ClearSubstitutions <: Function end

function (f::ClearSubstitutions)(output::Vector{Substitution}, readseq::LongSequence{DNAAlphabet{2}})
    return empty!(output)
end

"""
A simple function object that can be used with [`edit_substitutions`](@ref) or
[`edit_substitutions!`](@ref).

Randomly applies sequencing errors - substitutions - to all reads according to a
fixed, uniform, per-base probability.
"""
struct FixedProbSubstitutions <: Function
    prob::Float64
end
function (f::FixedProbSubstitutions)(output::Vector{Substitution}, readseq::LongSequence{DNAAlphabet{2}})
    p = f.prob
    dice = rand(length(readseq))
    dosub = dice .< p
    @inbounds for pos in 1:length(readseq)
        if dosub[pos]
            basenuc = reinterpret(UInt8, readseq[pos])
            push!(output, Substitution(pos, rand(posbases[basenuc])))
        end
    end
end

"""
    edit_substitutions!(f::Function, reads::Reads)

Add or remove sequencing errors to reads in the form of single base substitutions,
according to so some error generating or deleting function `f`.

!!! note
    This is the in-place, mutating version of [`edit_substitutions`](@ref), and
    since Pseudoseq's experted sequencing API is functional, this version is not
    exported. 

The function provided must be a function that accepts two arguments:
    1. A vector of Substitution structs, which the function will mutate.
    2. The nucleotide sequence of the read.

The function needs not return anything, and anything it does return will not be
used.

Some acceptable functions already provided with Pseudoseq include [`FixedProbSubstitutions`](@ref)
and [`ClearSubstitutions`](@ref). If you want to develop your own function, their
implementations provide a simple example of how it can be done. You can of course
also just create an anonymous function.
"""
function edit_substitutions!(f::Function, reads::Reads)
    vs = views(reads)
    subs = substitutions(reads)
    subsbuf = Vector{Substitution}()
    for (i, v) in enumerate(vs)
        readseq = extract_sequence(genome(reads), v)
        readsubs = get(subs, i, subsbuf)
        f(readsubs, readseq)
        if isempty(readsubs)           # Read no longer has any subs.
            delete!(subs, i)
        elseif readsubs === subsbuf    # Read did not have subs but now does.
            subs[i] = copy(subsbuf)
            empty!(subsbuf)
        else                           # Read did have subs, and still does.
            subs[i] = readsubs
        end
    end
end

"""
    edit_substitutions(f::Function, reads::Reads)

Add or remove sequencing errors to reads in the form of single base substitutions,
according to so some error generating or deleting function `f`.

The function provided must be a function that accepts two arguments:
    1. A vector of Substitution structs, which the function will mutate.
    2. The nucleotide sequence of the read.

The function needs not return anything, and anything it does return will not be
used.

Some acceptable functions already provided with Pseudoseq include [`FixedProbSubstitutions`](@ref)
and [`ClearSubstitutions`](@ref). If you want to develop your own function, their
implementations provide a simple example of how it can be done. You can of course
also just create an anonymous function.
"""
function edit_substitutions(f::Function, reads::Reads)
    newreads = typeof(reads)(genome(reads), views(reads), copy(substitutions(reads)))
    edit_substitutions!(f, newreads)
    return newreads
end

# Read file generation
# --------------------

@inline function _add_subs!(seq::BioSequence, subs::Vector{Substitution})
    for sub in subs
        seq[sub.pos] = sub.base
    end
end

function generate(wtr::FASTQ.Writer, reads::Reads{UnPaired,BasicSequencingView})
    vs = views(reads)
    g = genome(reads)
    subs = substitutions(reads)
    n = 0
    for i in eachindex(vs)
        v = vs[i]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(g, v)
        if haskey(subs, i)
            _add_subs!(seq, subs[i])
        end
        # Make the quality string for the read, (currently no real meaning!).
        qual = fill(Int8(30)::Int8, length(seq))
        #qual = Vector{Int}(undef, length(seq))
        #FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        # Create the name for the read...
        fragname = string("Refseq_", seqid(v))
        fqread = FASTQ.Record(fragname, seq, qual)
        write(wtr, fqread)
        n += 1
    end
    @info string("- ✔ Wrote ", n, " single end reads to FASTQ file")
    return reads
end

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{Paired,BasicSequencingView})
    vs = views(reads)
    g = genome(reads)
    subs = substitutions(reads)
    n = 0
    for i in eachindex(vs)
        v = vs[i]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(g, v)
        if haskey(subs, i)
            _add_subs!(seq, subs[i])
        end
        # Make the quality string for the read, (currently no real meaning!).
        #qual = Vector{Int}(undef, length(seq))
        #FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        qual = fill(Int8(30)::Int8, length(seq))
        # Create the name for the read...
        if isodd(i)
            fragname = string("readpair_", i)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R1W, fqread)
        else
            fragname = string("readpair_", i - 1)
            fqread = FASTQ.Record(fragname, seq, qual)
            write(R2W, fqread)
        end
        n += 1
    end
    @info string("- ✔ Wrote ", n, " paired end reads to FASTQ file")
    return reads
end

function prepare_tags(reads::Reads{Paired,TaggedSequencingView})
    tags = keys(summarize_tags(views(reads)))
    tagrange = 0x0000000000000000:0x00000000ffffffff
    tagseqs = DNAKmer{16}.(sample_values(tagrange, length(tags)))
    tagdict = Dict(zip(tags, tagseqs))
    return tagdict
end

function generate(R1W::FASTQ.Writer, R2W::FASTQ.Writer, reads::Reads{Paired,TaggedSequencingView})
    vs = views(reads)
    g = genome(reads)
    subs = substitutions(reads)
    
    bufseq = LongSequence{DNAAlphabet{2}}(rand(ACGT, 7))
    tagdict = prepare_tags(reads)
    n = 0
    
    for i in eachindex(vs)
        v = vs[i]
        # Extract subsequence from the reference and add errors.
        seq = extract_sequence(g, v)
        if haskey(subs, i)
            _add_subs!(seq, subs[i])
        end
        
        if isodd(i) # We are dealing with an R1 read
            tagseq = LongSequence{DNAAlphabet{2}}(tagdict[tag(v)])
            seq = tagseq * bufseq * seq
        end
        
        # Make the quality string for the read, (currently no real meaning!).
        qual = fill(Int8(30)::Int8, length(seq))
        #qual = Vector{Int}(undef, length(seq))
        #FASTQ.encode_quality_string!(FASTQ.SANGER_QUAL_ENCODING, fill(30, length(seq)), qual, 1, length(seq))
        
        # Create the name for the read...
        pairnum = Int(ceil(i / 2))
        fragname = string("readpair_", pairnum)
        fqread = FASTQ.Record(fragname, seq, qual)
        current_writer = ifelse(isodd(i), R1W, R2W)
        write(current_writer, fqread)
        
        n += 1
    end
    @info string("- ✔ Wrote ", n, " tagged paired end reads to FASTQ file")
    return reads
end

"""
    generate(R1name::String, R2name::String, reads::Reads{Paired,<:AbstractSequencingView})
    
This method only works for paired reads. Instead of interleaving R1 and R2 reads
in a single FASTQ file, R1 and R2 reads are partitioned into two seperate FASTQ
files.
"""
function generate(R1name::String, R2name::String, reads::Reads{Paired,<:AbstractSequencingView})
    R1W = FASTQ.Writer(open(R1name, "w"), quality_header = false)
    R2W = FASTQ.Writer(open(R2name, "w"), quality_header = false)
    try
        generate(R1W, R2W, reads)
    finally
        close(R1W)
        close(R2W)
    end
    return reads
end

generate(wtr::FASTQ.Writer, reads::Reads{Paired,<:AbstractSequencingView}) = generate(wtr, wtr, reads)

"""
    generate(filename::String, reads::Reads)
    
Write the `reads` out to a FASTQ formatted file with the given `filename`.

If this method is used with a paired-end read type, then the FASTQ file will
be interleaved; all R1 reads will be odd records, and all R2 reads will be even
records in the file.

!!! note
    Reads are named according to the sequence in the input genome they came
    from. e.g. `@Reference_1_R1` means the first sequence in the genome, and
    `@Reference_2_R1` means the second sequence in the genome.
"""
function generate(filename::String, reads::Reads)
    FASTQ.Writer(open(filename, "w"), quality_header = false) do wtr
        generate(wtr, reads)
    end
    return reads
end

function substitution_hist(reads::Reads{Paired,BasicSequencingView})
    hist = zeros(summarize_lengths(views(reads))[1])
    for i in Base.OneTo(nreads(reads))
        for s in get(substitutions(reads), i, Substitution[])
            hist[s.pos] += 1
        end
    end
    return hist
end