```@meta
# Information for Documenter
CurrentModule = GCATBase
```

```@contents
Pages = ["index.md"]
```

# Tutorial

## Pre-defined tuples

First, we need to load the required packages:

```@example rt
using GCATBase, BioSequences
```

We can access predefined sets of k-mers like codons, dinucleotides, and tetranucleotides:

```@example rt
GCATBase.codons
```

```@example rt
GCATBase.dinucs
```

```@example rt
GCATBase.tetranucs
```

## Tuples from sequences

A sequence can be split into non-overlapping tuples of a given length using the
[`GCATBase.split`](@ref) function. Here, we generate a random DNA sequence of
length 10 and split it into codons (tuples of length 3):

```@example rt
seq = randseq(DNAAlphabet{4}(), 10)
println(seq)
codons = split(seq; l=3)
println(codons)
```

## Shifting sequences

```@example rt
circshift(dna"AGCT")
```

Let's shift some RNA in the other direction:

```@example rt
circshift(rna"AGCU"; k=-1)
```

We can also use `circshift` together with `tuples` to generate all
cyclic permutations of k-mers in a sequence:

```@example rt
seq = rna"CAGCUUGAG"
join(circshift.(split(seq, l=3)))
```

# API

```@autodocs
Modules = [GCATBase]
Order   = [:type, :function]
```
