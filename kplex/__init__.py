"""k-plexity: the fraction of unique k-mers vs k, straight from an assembly.

Three commands:
  kplex run   FASTA -> k-plexity curve CSV (via FASTK/Kplex or KMC)
  kplex fit   curve CSV -> double-sigmoid parameters (L1, L2, s1, s2, k0_1, k0_2)
  kplex plot  curve CSV (+ optional fit) -> figure
"""
__version__ = "0.1.0"
