#!/usr/bin/env Rscript

### title: Print out heatmaps of IG combination usage in B cells author: Yiping
### Wang date: 02/08/2021

groups = c("cov", "ctr")
for (z in 1:length(groups)) {
    data_lungs_all_bcells = readRDS("/data/data_lungs_all_bcells.rds")
    data_lungs_all_bcells = subset(data_lungs_all_bcells, group == groups[z])
    uniqueIDs = unique(data_lungs_all_bcells$ID)
    
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
        "#D55E00", "#CC79A7")
    
    # read in data into IGtable variable, filtering out rows with NAs for both heavy
    # and light chains
    IGtable = read.table(paste0("table_IG_new_lungs_all_b_cell_int_all_chains_", 
        groups[z], ".csv"), header = T, sep = ",", quote = "")
    IGtable = IGtable[!is.na(IGtable$light) & !is.na(IGtable$heavy), ]
    
    # read in order of heavy and light chains from IG combination frequency heatmap
    uniquelight = unique(IGtable$light)
    uniqueheavy = unique(IGtable$heavy)
    lightorder = rev(read.table("lightorder.txt", header = F, quote = "")$V1)
    heavyorder = read.table("heavyorder.txt", header = F, quote = "")$V1
    uniquelight = lightorder
    uniqueheavy = heavyorder
    
    # store indices of where IGtable light and heavy chains match uniquelight and
    # uniqueheavy
    IGtable$lightidxs = match(IGtable$light, uniquelight)
    IGtable$heavyidxs = match(IGtable$heavy, uniqueheavy)
    
    # make list of all light and heavy chain combinations that occur in IGtable,
    # store in combs make list of combinations that are occur more than once, store
    # in dupcombs
    combs = c()
    dupcombs = c()
    for (i in 1:length(IGtable$light)) {
        if (sum(combs == paste0(IGtable$light[i], " ", IGtable$heavy[i])) != 0 & 
            sum(dupcombs == paste0(IGtable$light[i], " ", IGtable$heavy[i])) == 0) {
            dupcombs = append(dupcombs, paste0(IGtable$light[i], " ", IGtable$heavy[i]))
        }
        combs = append(combs, paste0(IGtable$light[i], " ", IGtable$heavy[i]))
    }
    
    # for duplicated combinations, it is possible that each entry may have different
    # values in the constant field find all unique values in constant field for each
    # duplicated combination in the cell for each duplicated combination, subdivide
    # the cell along the light chain axis, so that each subcell will represent a
    # unique constant value store the endpoints of each subcell in lightidxs and
    # lightidxsceil fields
    duptable = data.frame(light = character(), heavy = character(), constant = character(), 
        lightidxs = integer(), lightidxsceil = integer(), heavyidxs = integer())
    for (i in 1:length(dupcombs)) {
        words = strsplit(dupcombs[i], " ")[[1]]
        smalltable = na.omit(IGtable[IGtable$light == words[1] & IGtable$heavy == 
            words[2], ])
        smalltableconstant = unique(smalltable$constant)
        print(length(smalltableconstant))
        if (length(smalltableconstant) > 1) {
            for (j in 1:length(smalltableconstant)) {
                duptable = rbind(duptable, data.frame(light = words[1], heavy = words[2], 
                  constant = smalltableconstant[j], lightidxs = match(words[1], uniquelight) + 
                    (j - 1)/length(smalltableconstant), lightidxsceil = match(words[1], 
                    uniquelight) + j/length(smalltableconstant), heavyidxs = match(words[2], 
                    uniqueheavy)))
            }
        }
    }
    
    # create heatmap of all light ahd heavy chain combinations, with subcells for
    # multiple constant values appearing in one combination Figure 2g if using covid
    # samples Extended Data Figure B-cells f if using control samples
    p = ggplot(IGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = constant, 
        width = 1, height = 1)) + geom_tile()
    if (dim(duptable)[1] != 0) {
        p = p + geom_rect(data = duptable, aes(xmin = heavyidxs, xmax = heavyidxs + 
            1, ymin = lightidxs, ymax = lightidxsceil, fill = constant), colour = "black", 
            size = 0.3)
    }
    p = p + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 1), 
        labels = append(uniqueheavy, "")) + scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 
        1), labels = append(uniquelight, "")) + scale_fill_manual(breaks = c("IGHA1", 
        "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM"), values = cbPalette[1:8]) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
            axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.grid = element_line(size = 0.1, 
                colour = "black"), panel.ontop = T)
    ggsave(paste0("IG_combination_unique_", groups[z], ".pdf"), height = 12)
    
    # count number of times each IG chain combination occurs in IGtable, store in
    # sharedIGtable, along with lightidxs and heavyidxs
    sharedIGtable = data.frame(heavy = character(), light = character(), occurrences = integer(), 
        heavyidxs = integer(), lightidxs = integer())
    for (i in 1:length(uniqueIDs)) {
        IGtablesmall = IGtable[IGtable$barcode %in% names(data_lungs_all_bcells$ID[data_lungs_all_bcells$ID == 
            uniqueIDs[i]]), ]
        if (dim(IGtablesmall)[1] != 0) {
            checktable = data.frame(heavy = character(), light = character())
            for (j in 1:dim(IGtablesmall)[1]) {
                if (!is.na(IGtablesmall$heavy[j]) & !is.na(IGtablesmall$light[j])) {
                  if (sum(checktable$heavy == IGtablesmall$heavy[j] & checktable$light == 
                    IGtablesmall$light[j]) == 0) {
                    smalltable = data.frame(heavy = IGtablesmall$heavy[j], light = IGtablesmall$light[j])
                    checktable = rbind(checktable, smalltable)
                    
                    matchingarr = (sharedIGtable$heavy == IGtablesmall$heavy[j] & 
                      sharedIGtable$light == IGtablesmall$light[j])
                    if (sum(matchingarr) == 0) {
                      smalltable = data.frame(heavy = IGtablesmall$heavy[j], light = IGtablesmall$light[j], 
                        occurrences = 1, heavyidxs = match(IGtablesmall$heavy[j], 
                          uniqueheavy), lightidxs = match(IGtablesmall$light[j], 
                          uniquelight))
                      sharedIGtable = rbind(sharedIGtable, smalltable)
                    } else {
                      sharedIGtable$occurrences[matchingarr] = sharedIGtable$occurrences[matchingarr] + 
                        1
                    }
                  }
                }
            }
        }
    }
    
    # print heatmap of number of times each IG chain combination occurs in IGtable
    # Figure 2f if using covid samples Extended Data Figure B-cells e if using
    # control samples
    sharedIGtable$occurrences = as.character(sharedIGtable$occurrences)
    ggplot(sharedIGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = occurrences, 
        width = 1, height = 1)) + geom_tile() + scale_x_continuous("Heavy Chain", 
        breaks = 1:(length(uniqueheavy) + 1), labels = append(uniqueheavy, "")) + 
        scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 1), labels = append(uniquelight, 
            "")) + scale_fill_manual(breaks = c("1", "2", "3", "4", "5"), values = cbPalette[1:5]) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
            axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.grid = element_line(size = 0.1, 
                colour = "black"), panel.ontop = T)
    ggsave(paste0("IG_combination_unique_shared_across_samples_", groups[z], ".pdf"), 
        height = 12)
    
    # create new dataframe that stores information on IG combinations across covid
    # and control combinations are labeled as either found in COVID-19 or Control alone, or
    # shared among the two
    if (z == 1) {
        covandctr_sharedIGtable = sharedIGtable
        covandctr_sharedIGtable$occurrences = NULL
        covandctr_sharedIGtable$group = groups[z]
    } else {
        for (i in 1:dim(sharedIGtable)[1]) {
            if (sum(sharedIGtable$heavy[i] == covandctr_sharedIGtable$heavy & sharedIGtable$light[i] == 
                covandctr_sharedIGtable$light) != 0) {
                covandctr_sharedIGtable$group[sharedIGtable$heavy[i] == covandctr_sharedIGtable$heavy & 
                  sharedIGtable$light[i] == covandctr_sharedIGtable$light] = "shared"
            } else {
                smalltable = data.frame(heavy = sharedIGtable$heavy[i], light = sharedIGtable$light[i], 
                  heavyidxs = match(sharedIGtable$heavy[i], uniqueheavy), lightidxs = match(sharedIGtable$light[i], 
                    uniquelight), group = groups[z])
                covandctr_sharedIGtable = rbind(covandctr_sharedIGtable, smalltable)
            }
        }
    }
    
    # create heatmaps of IG combinations that occur in each individual sample
    for (i in 1:length(uniqueIDs)) {
        if (groups[z] == "ctr") {
            suffix = "ctr"
        } else {
            suffix = "cov"
        }
        IGtablesmall = IGtable[IGtable$barcode %in% names(data_lungs_all_bcells$ID[data_lungs_all_bcells$ID == 
            uniqueIDs[i]]), ]
        
        if (dim(IGtablesmall)[1] != 0) {
            combssmall = c()
            dupcombssmall = c()
            duptablesmall = data.frame(light = character(), heavy = character(), 
                constant = character(), lightidxs = integer(), lightidxsceil = integer(), 
                heavyidxs = integer())
            for (j in 1:length(IGtablesmall$light)) {
                if (sum(combssmall == paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j])) != 
                  0 & sum(dupcombssmall == paste0(IGtablesmall$light[j], " ", IGtablesmall$heavy[j])) == 
                  0) {
                  dupcombssmall = append(dupcombssmall, paste0(IGtablesmall$light[j], 
                    " ", IGtablesmall$heavy[j]))
                }
                combssmall = append(combssmall, paste0(IGtablesmall$light[j], " ", 
                  IGtablesmall$heavy[j]))
            }
            if (length(dupcombssmall) != 0) {
                for (k in 1:length(dupcombssmall)) {
                  words = strsplit(dupcombssmall[k], " ")[[1]]
                  smalltable = na.omit(IGtablesmall[IGtablesmall$light == words[1] & 
                    IGtablesmall$heavy == words[2], ])
                  smalltableconstant = unique(smalltable$constant)
                  if (length(smalltableconstant) > 1) {
                    for (l in 1:length(smalltableconstant)) {
                      duptablesmall = rbind(duptablesmall, data.frame(light = words[1], 
                        heavy = words[2], constant = smalltableconstant[l], lightidxs = match(words[1], 
                          uniquelight) + (l - 1)/length(smalltableconstant), lightidxsceil = match(words[1], 
                          uniquelight) + l/length(smalltableconstant), heavyidxs = match(words[2], 
                          uniqueheavy)))
                    }
                  }
                }
            }
            
            p = ggplot(IGtablesmall, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, 
                fill = constant, width = 1, height = 1)) + geom_tile()
            if (dim(duptablesmall)[1] > 0) {
                p = p + geom_rect(data = duptablesmall, aes(xmin = heavyidxs, xmax = heavyidxs + 
                  1, ymin = lightidxs, ymax = lightidxsceil, fill = constant), colour = "black", 
                  size = 0.3)
            }
            p = p + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 
                1), labels = append(uniqueheavy, ""), limits = c(1, length(uniqueheavy))) + 
                scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 
                  1), labels = append(uniquelight, ""), limits = c(1, length(uniquelight))) + 
                scale_fill_manual(breaks = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", 
                  "IGHG3", "IGHG4", "IGHGP", "IGHM"), values = cbPalette[1:8]) + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(r = 0)), 
                  axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), 
                  panel.grid.minor = element_blank(), panel.background = element_blank(), 
                  panel.grid = element_line(size = 0.1, colour = "black"), panel.ontop = T)  # + xlim(0,100) + ylim(0,100)
            ggsave(paste0("IG_combination_unique_", uniqueIDs[i], "_", suffix, ".pdf"), 
                height = 12)
            print("THERE")
        }
    }
    # scale_fill_manual(breaks = c('IGHA1', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4',
    # 'IGHM'),values=cbPalette[1:6]) +
}

# plot IG combinations that are found in COVID-19 alone, Control alone, or in both
ggplot(covandctr_sharedIGtable, aes(x = heavyidxs + 0.5, y = lightidxs + 0.5, fill = group, 
    width = 1, height = 1)) + geom_tile() + scale_x_continuous("Heavy Chain", breaks = 1:(length(uniqueheavy) + 
    1), labels = append(uniqueheavy, "")) + scale_y_continuous("Light Chain", breaks = 1:(length(uniquelight) + 
    1), labels = append(uniquelight, "")) + scale_fill_manual(breaks = c("cov", "ctr", 
    "shared"), values = cbPalette[1:3]) + theme(axis.text.x = element_text(angle = 90, 
    hjust = 1, margin = margin(r = 0)), axis.text.y = element_text(vjust = 0, margin = margin(r = 0)), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), panel.grid = element_line(size = 0.1, 
        colour = "black"), panel.ontop = T)
ggsave(paste0("Extended Data Figure B-cells g.pdf"), height = 12)
