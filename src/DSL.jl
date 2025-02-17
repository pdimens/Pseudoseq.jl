# DSL
# ===

# Types and methods that make up the Pseusoseq DSL API.

struct Amplifier{F<:Function}
    pred::F
    n::Int
end
amplify(n::Int) = Amplifier(x -> true, n)
amplify(pref::Function, n::Int) = Amplifier(pred, n)
(a::Amplifier)(p::Molecules) = amplify(a.pred, p, a.n)
Base.show(io::IO, amp::Amplifier) = print(io, "Amplifier(x", amp.n, ")")

struct Fragmenter
    meansize::Int
end
Base.show(io::IO, frag::Fragmenter) = print(io, "Fragmenter(mean size: ", frag.meansize, "bp)")

"""
    fragment(meansize::SequenceLength)

Construct a Fragmenter transform that fragments any `Molecules` passed to it, to
an average length of `meansize`.
"""
fragment(meansize::SequenceLength) = Fragmenter(meansize.val)
(f::Fragmenter)(p::Molecules) = fragment(p, f.meansize)

struct Tagger
    ntags::Int
end
Base.show(io::IO, tag::Tagger) = print(io, "Tagger(", tag.ntags, " molecular tags)")

tag(ntags::Int) = Tagger(ntags)
(t::Tagger)(p::Molecules) = tag(p, t.ntags)

struct NSubSampler
    nsamples::Int
end
subsample(n::Int) = NSubSampler(n)
(s::NSubSampler)(p::Molecules) = subsample(p, s.nsamples)

struct CovSubSampler{T<:Union{SequenceLength,Tuple{SequenceLength,SequenceLength}}}
    cov::ExpectedCoverage
    readlen::T
end
Base.show(io::IO, covsampler::CovSubSampler) = print(io, "Coverage Subsampler\n  coverage: ", covsampler.cov.val, "X\n  length: 2 * ", covsampler.readlen[1])

subsample(cov::ExpectedCoverage, rlens)  = CovSubSampler{typeof(rlens)}(cov, rlens)
(s::CovSubSampler{T})(p::Molecules) where {T} = subsample(p, s.cov, s.readlen)

struct Selector{F<:Function}
    f::F
end
Base.show(io::IO, sel::Selector) = print(io, "Selector(", sel.f, ")")

select(f::Function) = Selector(f)
(s::Selector{F})(p::Molecules) where {F<:Function} = select(s.f, p)

makereads(len::SequenceLength) = unpaired_reads(len.val)
makereads() = unpaired_reads(nothing)
makereads(lena::SequenceLength, lenb::SequenceLength) = paired_reads(lena.val, lenb.val)
makereads(lens::Tuple{SequenceLength,SequenceLength}) = makereads(lens...)


struct PairedReads
    flen::Int
    rlen::Int
end
paired_reads(flen::Int, rlen::Int = flen) = PairedReads(flen, rlen)
(pr::PairedReads)(p::Molecules) = paired_reads(p, pr.flen, pr.rlen)

struct UnPairedReads{T<:Union{Nothing,Int}}
    len::T
end
unpaired_reads(len::T) where {T<:Union{Nothing,Int}} = UnPairedReads(len)
(ur::UnPairedReads)(p::Molecules) = unpaired_reads(p, ur.len)

struct SubstitutionMaker{F<:Function}
    fun::F
end
Base.show(io::IO, submaker::SubstitutionMaker) = print(io, "SubstitutionMaker(", submaker.fun, ")")
make_substitutions(f::F) where {F<:Function} = SubstitutionMaker{F}(f)
(sm::SubstitutionMaker{F})(p::Reads) where {F<:Function} = edit_substitutions(sm.fun, p)

struct FileGenerator
    filename::String
end
generate(filename::String) = FileGenerator(filename)
(fg::FileGenerator)(r::Reads) = generate(fg.filename, r)

struct FilePairGenerator
    F1::String
    F2::String
end
generate(F1::String, F2::String) = FilePairGenerator(F1, F2)
(fg::FilePairGenerator)(r::Reads) = generate(fg.F1, fg.F2, r)