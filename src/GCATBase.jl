module GCATBase

export alltuples, tuples, circshift

using BioSequences, BioSymbols
using DocStringExtensions
using Pipe: @pipe

"""
    $(TYPEDSIGNATURES)

Create all possible tuples of length `l` from the given `alphabet`.
The tuples are returned as a vector of `LongDNA` or `LongRNA` sequences,
depending on the specified sequence type `S`.
"""
function alltuples(alphabet::Vector{T}, l::Int, S::Union{Type{DNA},Type{RNA}}) where {T<:BioSymbol}
    v = repeat([alphabet], l)
    ts = reshape(collect(Iterators.product(v...)), 1, :)[1, :] #  make vector
    if S == DNA
        return [LongDNA{4}(t) for t in ts]
    elseif S == RNA
        return [LongRNA{4}(t) for t in ts]
    else
        error("Unsupported sequence type. Must be DNA or RNA.")
    end
end

const dinucs = alltuples([DNA_A, DNA_T, DNA_C, DNA_G], 2, DNA)
const codons = alltuples([DNA_A, DNA_T, DNA_C, DNA_G], 3, DNA)
const tetranucs = alltuples([DNA_A, DNA_T, DNA_C, DNA_G], 4, DNA)

import Base.split

"""
    $(TYPEDSIGNATURES)

Create tuples of length `l` by splitting a sequence `seq`.
`seq` must implement the `collect` method.
"""
function Base.split(seq::T; l=3) where {T<:BioSequence}
    S = typeof(seq)
    s = collect(seq) # make vector
    n = length(s) # length of sequence
    e = (n ÷ l) * l # end of sequence as a multiple of l
    @pipe s[1:e] |> reshape(_, l, :) |> eachcol |> map(S, _)
end

import Base.circshift

"""
    $(TYPEDSIGNATURES)

Shift the elements in a sequence for one position to the left (AUGC -> UGCA).
"""
function Base.circshift(seq::T; k=1) where {T<:BioSequence}
    l = length(seq)
    h = k % l
    i = h > 0 ? h : l + h
    if i > 0 # Anything to do at all?
        S = typeof(seq)
        a = collect(seq)
        S([a[(i+1):end]..., a[1:i]...])
    else
        seq
    end
end

end # module end