module Bioinfomatics where

data Nucleotide = A | C | G | T
                  deriving (Eq, Show, Read)

-- Gets the complement of a Nucleotide.
compn :: Nucleotide -> Nucleotide
compn A = T
compn C = G
compn G = C
compn T = A

type Sequence = [Nucleotide]

-- Reads a Sequence from a String.
readSequence :: String -> Sequence
readSequence ss = map (read . (:[])) ss :: [Nucleotide]

-- Gets the complement of a Sequence.
comps :: Sequence -> Sequence
comps = map compn

-- Gets the reverse complement of a Sequence.
revcomps :: Sequence -> Sequence
revcomps = reverse . comps

data Codon = Codon Nucleotide Nucleotide Nucleotide

data AminoAcid = Ala | Arg | Asn | Asp | Cys | Gln | Glu | Gly | His | Ile |
                 Leu | Lys | Met | Phe | Pro | Ser | Thr | Trp | Tyr | Val |
                 STOP

-- Translates a Codon into an AminoAcid.
translatec :: Codon -> AminoAcid
translatec (Codon A A A) = Lys
translatec (Codon A A C) = Asn
translatec (Codon A A G) = Lys
translatec (Codon A A T) = Asn

translatec (Codon A C A) = Thr
translatec (Codon A C C) = Thr
translatec (Codon A C G) = Thr
translatec (Codon A C T) = Thr

translatec (Codon A G A) = Arg
translatec (Codon A G C) = Ser
translatec (Codon A G G) = Arg
translatec (Codon A G T) = Ser

translatec (Codon A T A) = Ile
translatec (Codon A T C) = Ile
translatec (Codon A T G) = Met
translatec (Codon A T T) = Ile

translatec (Codon C A A) = Gln
translatec (Codon C A C) = His
translatec (Codon C A G) = Gln
translatec (Codon C A T) = His

translatec (Codon C C A) = Pro
translatec (Codon C C C) = Pro
translatec (Codon C C G) = Pro
translatec (Codon C C T) = Pro

translatec (Codon C G A) = Arg
translatec (Codon C G C) = Arg
translatec (Codon C G G) = Arg
translatec (Codon C G T) = Arg

translatec (Codon C T A) = Leu
translatec (Codon C T C) = Leu
translatec (Codon C T G) = Leu
translatec (Codon C T T) = Leu

translatec (Codon G A A) = Glu
translatec (Codon G A C) = Asp
translatec (Codon G A G) = Glu
translatec (Codon G A T) = Asp

translatec (Codon G C A) = Ala
translatec (Codon G C C) = Ala
translatec (Codon G C G) = Ala
translatec (Codon G C T) = Ala

translatec (Codon G G A) = Gly
translatec (Codon G G C) = Gly
translatec (Codon G G G) = Gly
translatec (Codon G G T) = Gly

translatec (Codon G T A) = Val
translatec (Codon G T C) = Val
translatec (Codon G T G) = Val
translatec (Codon G T T) = Val

translatec (Codon T A A) = STOP
translatec (Codon T A C) = Tyr
translatec (Codon T A G) = STOP
translatec (Codon T A T) = Tyr

translatec (Codon T C A) = Ser
translatec (Codon T C C) = Ser
translatec (Codon T C G) = Ser
translatec (Codon T C T) = Ser

translatec (Codon T G A) = STOP
translatec (Codon T G C) = Cys
translatec (Codon T G G) = Trp
translatec (Codon T G T) = Cys

translatec (Codon T T A) = Leu
translatec (Codon T T C) = Phe
translatec (Codon T T G) = Leu
translatec (Codon T T T) = Phe

-- Gets a list of Codons from a Sequence.
codons :: Sequence -> [Codon]
codons (n1:n2:n3:ns) = Codon n1 n2 n3 : codons ns
codons _ = []

-- Gets the dotplot of two sequences.
dotplot :: Sequence -> Sequence -> [[Bool]]
dotplot (x:xs) ys = map (==x) ys : dotplot xs ys
dotplot _      _  = []

-- Translates a Sequence into a list of AminoAcid, per an integer frame.
-- Frame should be 0, 1, or 2, but no such restriction is necessary here.
type Frame = Int
translates :: Frame -> Sequence -> [AminoAcid]
translates n = map translatec . drop n . codons
