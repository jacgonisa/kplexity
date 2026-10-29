"""k-plexity: the fraction of unique k-mers vs k, straight from an assembly.

Three commands:
  kplex run   FASTA -> k-plexity curve CSV (via FASTK/Kplex or KMC)
  kplex fit   curve CSV -> model parameters (A, S, lam/divergence, k0, b; or the double sigmoid)
  kplex plot  curve CSV (+ optional fit) -> figure
"""
__version__ = "0.3.0"
