  source("OPAM.library.v7.R")     # ssGSEA library

   OPAM.project.dataset.5(
      input.ds                 = "~/CGP2012/CCLE/CCLE_Expression_Entrez_2012-04-06_plus_Achilles_appendix.gct",
      output.ds                = "~/CGP2013/signatures/RAS/CCLE_MSigDB.PATHWAYS.gct",
      gene.set.databases       =  c("~/CGP2013/signatures/RAS/c2.all.v3.1.symbols.gmt",
                                    "~/CGP2013/signatures/RAS/c3.all.v3.1.symbols.gmt",
                                    "~/CGP2013/signatures/RAS/c5.all.v3.1.symbols.gmt",
                                    "~/CGP2013/signatures/RAS/c6.all.v3.1.symbols.gmt"),
      gene.set.selection       = "ALL",
      sample.norm.type         = "rank",  # "rank", "log" or "log.rank"
      weight                   = 0.75,
      statistic                = "area.under.RES",
      output.score.type        = "ES",
      combine.mode             = "combine.add",  # "combine.off", "combine.replace", "combine.add"
      nperm                    =  1,
      min.overlap              =  1,
      correl.type              =  "z.score")             # "rank", "z.score", "symm.rank"