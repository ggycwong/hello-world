#Code for:
#Figure 8A. The number of A. thaliana specific, Brassicaceae specific, and conserved miRNA families identified with and without Major targets
#Figure 8B. The number of Major targets identififed for A. thaliana specific, Brassicaceae specific, and conserved miRNA families

library(dplyr)
library(magicfor)
library(ggplot2)


#The 'miRNAome' all miRNAs retrieved from miRbase v22 used for analysis
mirna <- read.csv("Chp3_Results/Figure_8_9/Arabidopsis_miRNome_categorised.csv")
allmirna <- unique(mirna$miRNA) #All A. thaliana miRNAs annotated on miRBase v22

#psRNATarget predicted targets of all miRNAs analysed
#Note: takes longer to load
psRNA_data <- read.delim("Chp3_Results/Figure_8_9/psRNATarget_ALL_miRBase_AthMiRNAs.txt")

#Degradome sequencing data for all miRNAs analysed
#Note: takes longer to load
deg_raw <- read.csv("Chp3_Results/Figure_8_9/ALL_miRBase_AthMiRNAs_Verified_Targets_raw.csv")
deg_data <- deg_raw

lib_co <- 10 #A library cut-off of 10% is used. Library cut-off can be changed manually

magic_for(func = put, silent = TRUE)

for(i in 1:length(allmirna)) { #Loops through all miRNAs

  
  print("########################################")
  print(paste("Process", i, allmirna[i], sep = ":"))
  
  mirna_data <- subset(deg_data, deg_data$Small_RNA %in% allmirna[i])
  
  ###Parameter B and C: The Cleavage tag abundance and miRNA target Category parameter
  #Validated targets of miR162, miR396, miR398 and miR408 results in a 1nt bulge in the target which means it can only be detected by the degradome pipeline if a 1nt offset is included
  if(any(grepl(paste(c("miR162", "miR396", "miR398", "miR408"),collapse="|"), mirna_data$Small_RNA)) #Factor in a 1nt offset for miR162, miR396, miR398 and miR408
  ) {
    shift_data <- mirna_data 
  } else {
    shift_data <- subset(mirna_data, mirna_data$Shift == "0") #Ignores the 1nt offset for all other miRNAs
  }
  
  ###Parameter B and C: The Cleavage tag abundance and miRNA target Category parameter
  ord_deg <- shift_data %>%
    subset(Category != "Cat_3" & Category != "Cat_4") %>%  #Exclude Category 3 and 4 target data
    subset(Deg_count >=5) %>%  #Include cleavage tags that are >=5 
    arrange(-Deg_count)  #Order by cleavage tag abundance
  
  max_shift <- ord_deg %>% #This is for "miR396", "miR398", "miR408", "miR162". It takes the offset position with the highest cleavage tags so that only one offset position is taken into account
    group_by(Degradome, Target_gene) %>% 
    summarise(Deg_count = max(Deg_count) , na.rm = TRUE)
  
  count_data <- max_shift %>%
    group_by(Target_gene) %>%
    summarise(Library_count =n()) %>% #Count the number of libraries a gene occurs in
    arrange(-Library_count) %>% #Order by number of libraries
    as.data.frame()
  
  #Formatting transcript ID to gene ID
  count_data$Target_gene <- gsub( "\\..*", "", as.character(count_data$Target_gene)) #Make 'Target_gene' a gene ID instead of a transcript ID
  count_data <- count_data[!duplicated(count_data$Target_gene),]  #Remove duplicated target entries from transcript isomers
  
  ###Parameter D: Library Cut-off
  num_lib <- length(unique(deg_data$Degradome)) #34 A. thaliana degradome libraries are used for this analysis

  b_d_filtered <- subset(count_data,(count_data$Library_count >= (num_lib*(lib_co/100)))) #Targets filtered through parameters B-D
  print(b_d_filtered$Target_gene)
  
  ###Parameter A: The psRNATarget Expectation score cut-off
  ##Predicted targets----------------------------------------
  psRNA_data$Target_Acc. = gsub( "\\..*", "", as.character(psRNA_data$Target_Acc.)) #Make 'Target_Acc.' a gene ID instead of a transcript ID 
  psRNA_data <-  psRNA_data[!duplicated(psRNA_data$Target_Acc.),] #and get rid of duplicates
  
  if(any(grepl(paste(c("miR167",  "miR398",  "miR408"), collapse="|"), allmirna[i]))) {
    psRNA_exp <- subset(psRNA_data, psRNA_data$Expectation <= 5.0) #Expectation score cut-off of 5 if miR167, miR398, miR408.
  } else {
    psRNA_exp <- subset(psRNA_data, psRNA_data$Expectation <= 3.0) #Expectation score cut-off of 3 for all other miRNAs
  }
  
  #See results in the loop
  majtar_check <- as.data.frame(intersect(b_d_filtered$Target_gene, psRNA_exp$Target_Acc.))
  majtar_check$miRNA <- (replicate(nrow(majtar_check), allmirna[i]))
  names(majtar_check) <- c("ID", "miRNA") 
  print(majtar_check)
  #See results in the loop
  
  ID <- intersect(b_d_filtered$Target_gene, psRNA_exp$Target_Acc.) #These are the Major Targets
  if(length(ID) == 0) {
    ID <- print("No Major Targets") #If no Major targets are identified, print "No Major Target"
  }
  
  miRNA <- replicate(length(ID), allmirna[i]) #The miRNA
  Library_cutoff <- replicate(length(ID), lib_co) #The library cut-off
  miRNA_sequence <- replicate(length(ID), mirna_data$sRNA_seq[1]) #The miRNA sequence
  
  put(miRNA, ID, miRNA_sequence, Library_cutoff)
}
major_targets <- magic_result_as_vector() #Output for for loop results
majr_tar_df <- as.data.frame(major_targets)
head(majr_tar_df)

#Label miRNAs as guide or passenger strands
strands <- read.csv("Chp3_Results/Figure_8_9/miRNA_strands.csv")
label_strand <- left_join(majr_tar_df, strands)  #label guide/passenger strands.

#Label the conservation of the miRNA families
label_conservation <- left_join(label_strand, mirna) 
col_order <- c("miRNA", "miRNA_family", "ID", "Conservation", "miRNA_sequence", "Strand", "Library_cutoff") #Rearrange columns
labelled_mirna <- label_conservation[, col_order]

#write.csv(labelled_mirna, "Chp3_Results/Figure_8_9/Results_Output/OutputTable1.csv")


#-----------------------------------------------------------------------------------------------------

#IsomiRs in the same miRNA family will target many of the same targets 
#Remove all duplicated targets of the same miRNA family so that the target is only counted once per miRNA family
all_maj_targets <- labelled_mirna %>% 
  group_by(miRNA_family) %>% #Group by target miRNA family
  filter(!duplicated(ID)) #Remove duplicated targets for each miRNA family as they will have duplicated Major targets

#-----------------------------------------------------------------------------------------------------

#The number of miRNA families per conservation group
num_miRNA <- all_maj_targets %>%
  group_by(Conservation, Strand) %>%
  filter(!duplicated(miRNA_family)) %>% #Remove duplicated targets for each miRNA family as they will have duplicated Major targets
  count(Conservation, Strand, name = "No. miRNA Families") %>%
  print()

#-----------------------------------------------------------------------------------------------------

#The number of miRNA families with Major targets
num_miRNA_w_targets <- all_maj_targets %>%
  group_by(Conservation, Strand) %>%
  subset(ID != "No Major Targets") %>% #Only analyse miRNA with no Major targets
  filter(!duplicated(miRNA_family)) %>% #Remove duplicated miRNA sequences as isomiRs with the same sequence will have the same Major targets
  count(Conservation, Strand, name = "No. miRNA Families with Major Targets") %>%
  as.data.frame() %>%
  print()

#-----------------------------------------------------------------------------------------------------

#The number of Major targets per conservation group
num_targets_per_cons <- all_maj_targets %>%
  group_by(Conservation, Strand) %>%
  subset(ID != "No Major Targets") %>% #Remove miRNAs with no Major targets
  count(Conservation, name = "No. Major Targets") %>% #Count the number of Major targets per conservation group
  as.data.frame() %>%
  print()

#-----------------------------------------------------------------------------------------------------

#Figure 8b. Category score for each Major target

magic_for(func = put, silent = TRUE)

for(i in 1:nrow(all_maj_targets)) { #Find Category Score for all Major
  
  deg_data$Target_gene <- gsub( "\\..*", "", as.character(deg_data$Target_gene)) #Make 'Target_gene' a gene ID instead of a transcript ID
  print("###########################################################")
  print(paste(i, all_maj_targets[i,]$miRNA, all_maj_targets[i,]$ID, sep = ": "))
  
  mirna_df <- deg_data %>% 
    subset(x = ., all_maj_targets[i,]$ID == Target_gene) %>%
    subset(x = ., Deg_count >= 5) %>% #Exclude data with cleavage tag abundances equal to or less than 5
    subset(x = ., Category != "Cat_3" & Category != "Cat_4") %>% #Include only Category 1 and Category 2 targets
    arrange(Category, desc(Deg_count)) %>%
    filter(!duplicated(Degradome)) #We only want to count the library once
  
  ##Calculating the Category score
  ##This score takes into account the number of Category 1 and Category 2 libraries
  positive_libs <- nrow(mirna_df) #Number of libraries the Major target occurs in
  num_lib <- length(unique(deg_data$Degradome)) #34 A. thaliana degradome libraries are used for this analysis
  Cat_count <- count(mirna_df, Category)
  Cat_1 <- nrow(subset(mirna_df, mirna_df$Category == "Cat_1")) #Number of libraries where target is Category 1
  Cat_2 <- nrow(subset(mirna_df, mirna_df$Category == "Cat_2")) #Number of libraries where target is Category 2
  
  ##Equation is ('No. of libraries target is present in' / 'Total no. of libraries analysed') * ((('No. of Category 1 targets' * 5)  + ('No. of Category 2 targets'))/'No. of libraries target is present in')
  if(Cat_1 == 0){
    Category_score <- positive_libs / num_lib * ((Cat_2 * 1)/positive_libs) #Weighted score for when target is Category 1 and Category 2 (score of 3 and 1, respectively)
  }else{
    Category_score <- positive_libs / num_lib * (((Cat_1 * 5) + (Cat_2 * 1))/positive_libs) #Weighted score for when target is Category 1 and Category 2 (score of 3 and 1, respectively)
  }
  
  #Output columns for for loop.
  miRNA <- all_maj_targets[i,]$miRNA
  miRNA_family <- all_maj_targets[i,]$miRNA_family
  ID <- all_maj_targets[i,]$ID
  Conservation <- all_maj_targets[i,]$Conservation
  Category_score <- round(Category_score, digits = 3)
  Category <- mirna_df$Category[1]
  Cleavage_tag <- mirna_df$Deg_count[1]
  miRNA_sequence <- all_maj_targets[i,]$miRNA_sequence
  Strand <- all_maj_targets[i,]$Strand
  
  put(miRNA, miRNA_family, ID, Conservation, Category_score, Category, Cleavage_tag, miRNA_sequence, Strand)
}

#Targets with Category scores
cat <- magic_result_as_dataframe() #Output for for loop results
cat <- select(cat, -i)

#Supplementary Table 4. Major targets of miRNAs of all category groups with category scores
#write.csv(cat, row.names = FALSE, "Chp3_Results/Figure_8_9/Results_Output/SuppTable_4_test.csv")


#-----------------------------------------------------------------------------------------------------

#The number of Major targets with category 1 targets per conservation group
num_cat_1_targets_per_cons <- cat %>%
  subset(Category == "Cat_1") %>%
  group_by(Conservation, Strand) %>%
  subset(ID != "No Major Targets") %>% #Remove miRNAs with no Major targets
  count(Conservation, name = "No. Major Targets with Category 1 Targets") %>% #Count the number of Major targets per conservation group
  as.data.frame() %>%
  print()

#-----------------------------------------------------------------------------------------------------

#The average category score for each conservation group
cat_score_per_cons <- cat %>% 
  group_by(Conservation, Strand) %>% #Analyse per conservation group, and analyse between guide and passenger strands for conserved miRNAs
  summarize("Mean Category Score" = round(mean(Category_score, na.rm = TRUE),2)) %>% #Find mean category score and round to 2 decimal places)
  print()

#-----------------------------------------------------------------------------------------------------




#Figure 8A. The number of A. thaliana specific, Brassicaceae specific and conserved miRNA Families with Major targets

fig8a_df <- num_miRNA %>% 
  merge(num_miRNA_w_targets) %>% #Merge dataframe of No. of miRNA families with No. of miRNA families with Major targets
  mutate("No. miRNA Families without Major Targets" = (num_miRNA$`No. miRNA Families` - num_miRNA_w_targets[,3])) %>% #Also need No. of miRNA families without Major targets
  select(-"No. miRNA Families") #Remove No. miRNA Families to produce graph

fig8a_df$Conservation <- c("A. thaliana", "Brassicaceae", "Conserved (Guide)", "Conserved (Passenger)") #Change conservation group names

fig8a_df_long <- reshape2::melt(fig8a_df, value.name = "num_miRNAs") #Convert data.frame table to long form
print(fig8a_df_long)

#Graph 1
fig8a <- ggplot(data = fig8a_df_long, aes(x = Conservation, y = num_miRNAs, 
                                          fill = factor(variable, levels = c("No. miRNA Families without Major Targets", 
                                                                   "No. miRNA Families with Major Targets")))) +
  geom_bar(stat ="identity", position = "stack") +
  labs(x = "Conservation Group", y = "No. of miRNA Families") +
  theme(
    axis.text.y = element_text(size=10.5),
    axis.text.x = element_text(face = "bold", size = 8),
    axis.title=element_text(size=12,face="bold"),
    strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"),
    strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white')
  ) 

fig8a +
  scale_fill_manual(values=c("#CC6699", "#FF9966"),
                    name="miRNA Families",
                    breaks=c("No. miRNA Families with Major Targets", "No. miRNA Families without Major Targets"),
                    labels=c("with Major Targets", "without Major Targets")) +
  theme(legend.position="bottom")

#-----------------------------------------------------------------------------------------------------

#Figure 8B. The number and category score of Major targets for A. thaliana specific, Brassicaceae specific and conserved miRNA families

fig8b_df <- num_targets_per_cons #Number of Major targets for each conservation group

fig8b_df$Conservation <- c("A. thaliana", "Brassicaceae", "Conserved (Guide)", "Conserved (Passenger)") #Change conservation group names

fig8b_df_long <- reshape2::melt(fig8b_df, value.name = "No. Major Targets") #Convert data.frame table to long form
print(fig8b_df_long)

fig8b <- ggplot(data = fig8b_df_long, aes(x = Conservation, y = `No. Major Targets`)) +
  geom_bar(stat ="identity", fill = "#66CC99") +
  labs(x = "Conservation Group", y = "No. of Major Targets") +
  theme(
    axis.text.y = element_text(size=10.5),
    axis.text.x = element_text(face = "bold", size = 8),
    axis.title=element_text(size=12,face="bold"),
    strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"),
    strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white')
  ) 
print(fig8b)


#-----------------------------------------------------------------------------------------------------

#Figure 8C. The fold difference of predicted targets compared to Major targets for A. thaliana specific, Brassicaceae specific and conserved miRNA families

psRNA_cutoff <- subset(psRNA_data, psRNA_data$Expectation <= 3.0)
colnames(psRNA_cutoff)[1] <- "miRNA" #Change "miRNA_Acc." to "miRNA" to join with the "mirna" dataframe

#The total number of target predicted by psRNATarget for A. thaliana specific, Brassicaceae specific and conserved miRNA families
total_pred_targets <- psRNA_cutoff %>%
  left_join(mirna, by = "miRNA") %>% #label conservation group
  left_join(strands, by = "miRNA") %>% #label guide/passenger strand
  mutate(Target_Acc. = gsub("\\..*", "", Target_Acc.)) %>% #Change predicted targets to a gene ID instead of a transcript ID 
  filter(!duplicated(Target_Acc.)) %>% #Remove duplicated Major targets from isomers
  count(Conservation, Strand, name = "No. Predicted Targets") %>%
  print()

pred_maj_df <- merge(total_pred_targets, num_targets_per_cons) #Join number of predicted targets with number of Major targets
pred_maj_df$fold_diff <- round(pred_maj_df$`No. Predicted Targets`/pred_maj_df$`No. Major Targets`, 2) #Calculate fold different between number of predicted targets with number of Major targets
print(pred_maj_df)

pred_maj_df$Conservation <- c("A. thaliana", "Brassicaceae", "Conserved (Guide)", "Conserved (Passenger)") #Change conservation group names

fig8c <- ggplot(data = pred_maj_df, aes(x = Conservation, y = fold_diff)) +
  geom_bar(stat ="identity", fill = "#99CCFF") +
  labs(x = "Conservation Group", y = "No. of Predicted Targets/No. of Major Targets") +
  theme(
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(face = "bold", size = 8),
    axis.title=element_text(size=12,face="bold"),
    strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"),
    strip.text.x = element_text(size = 12, color = "white", face = "bold.italic"),
    panel.grid.major = element_line(colour = "grey85"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white')
  ) 
print(fig8c)



#-----------------------------------------------------------------------------------------------------

#Table 5. All the data above in one table - No. of miRNAs, No. of miRNAs with Major targets, No. of predicted targets, Number of category 1 targets, average Category score
cons_df <- num_miRNA %>% 
  merge(num_miRNA_w_targets) %>%
  merge(num_targets_per_cons) %>%
  merge(total_pred_targets) %>%
  merge(num_cat_1_targets_per_cons) %>%
  merge(cat_score_per_cons) %>%
  print()

cons_df$Conservation <- c("A. thaliana", "Brassicaceae", "Conserved (Guide)", "Conserved (Passenger)") #Change conservation group names


#-----------------------------------------------------------------------------------------------------
#Figure 9. Major targets identified for all A. thaliana miRNAs by conservation group


label_guide_pass <- mutate(labelled_mirna, Conservation = paste(labelled_mirna$Conservation, labelled_mirna$Strand, sep = "_")) #Label guide or passenger strand for conserved miRNAs

family_count <- label_guide_pass %>%
  group_by(miRNA_family) %>% #Group by target miRNA family and ID
  filter(!duplicated(ID)) %>% #Remove duplicated miRNA sequences as isomiRs with the same sequence will have the same Major targets
  subset(ID != "No Major Targets") %>%
  group_by(Conservation) %>%
  count(miRNA_family) %>%
  arrange(-n) %>%
  print()

data = family_count

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(as.factor(data$Conservation)), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$Conservation=rep(levels(as.factor(data$Conservation)), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(Conservation)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(Conservation) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
ggplot(data, aes(x=as.factor(id), y=n, fill=Conservation)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=n, fill=Conservation), stat="identity") +
  
  # Add a val=0/5/10/15 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 5, xend = start, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 15, xend = start, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the n of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0, 5, 10, 15), label = c("0", "5", "10", "15") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=n, fill=Conservation), stat="identity", alpha=0.5) +
  ylim(-20,25) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=n+2, label=miRNA_family, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -1, xend = end, yend = -1), colour = "black", alpha=1.5, size=1.0 , inherit.aes = FALSE )









#-----------------------------------------------------------------------------------------------------

#I'm trying to figure out if there is any merit in this analysis.....
#Actually, I think it does. We see fewer miRNAs with Major targets in A. thaliana specific between Figure 9A and B.
#This means many Major targets of the A. thaliana specific miRNAs had a Cat Max = 2 suggesting a lower confidence as miRNA targets.

#Figure 9b or 10. Major targets with at least on Category 1 target identified for all A. thaliana miRNAs by conservation group

cat1_max_targets <- cat %>%
  subset(Category == "Cat_1") %>%
  group_by(Conservation, Strand) %>%
  subset(ID != "No Major Targets") %>%
  as.data.frame()

label_guide_pass <- mutate(cat1_max_targets, Conservation = paste(cat1_max_targets$Conservation, cat1_max_targets$Strand, sep = "_")) #Label guide or passenger strand for conserved miRNAs

family_count <- label_guide_pass %>%
  group_by(miRNA_family) %>% #Group by target miRNA family and ID
  filter(!duplicated(ID)) %>% #Remove duplicated miRNA sequences as isomiRs with the same sequence will have the same Major targets
  subset(ID != "No Major Targets") %>%
  group_by(Conservation) %>%
  count(miRNA_family) %>%
  arrange(-n) %>%
  print()

data = family_count

# Set a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(as.factor(data$Conservation)), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$Conservation=rep(levels(as.factor(data$Conservation)), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(Conservation)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>% 
  group_by(Conservation) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
ggplot(data, aes(x=as.factor(id), y=n, fill=Conservation)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=n, fill=Conservation), stat="identity") +
  
  # Add a val=0/5/10/15 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 5, xend = start, yend = 5), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 15, xend = start, yend = 15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the n of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0, 5, 10, 15), label = c("0", "5", "10", "15") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=n, fill=Conservation), stat="identity", alpha=0.5) +
  ylim(-20,25) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=n+2, label=miRNA_family, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -1, xend = end, yend = -1), colour = "black", alpha=1.5, size=1.0 , inherit.aes = FALSE )









