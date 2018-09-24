# Set GO data
set_go_data <- compiler::cmpfun(function(organism, ontology) {
    species <- switch(organism, human = "org.Hs.eg.db",
                      fly = "org.Dm.eg.db",
                      mouse = "org.Mm.eg.db",
                      rat = "org.Rn.eg.db",
                      yeast = "org.Sc.sgd.db",
                      zebrafish = "org.Dr.eg.db",
                      worm = "org.Ce.eg.db",
                      arabidopsis = "org.At.tair.db",
                      ecolik12 = "org.EcK12.eg.db",
                      bovine = "org.Bt.eg.db",
                      canine = "org.Cf.eg.db",
                      anopheles = "org.Ag.eg.db",
                      ecsakai = "org.EcSakai.eg.db",
                      chicken = "org.Gg.eg.db",
                      chimp = "org.Pt.eg.db",
                      malaria = "org.Pf.plasmo.db",
                      rhesus = "org.Mmu.eg.db",
                      pig = "org.Ss.eg.db",
                      xenopus = "org.Xl.eg.db")
    GOSemSim::godata(species, ont = ontology, computeIC = TRUE)
})
