df = read.csv("~/Wang/output/Paper_figs/Fig5/Peaks_Gpr17_Cor_500KB_Stress_2024-02-12.csv")

endPost = as.numeric(as.character(gsub("^.*\\-","", df$peak )))

#startPos = gsub("\\chr-", df$peak )

endPos1 = 31951636 - endPost

df$Distance_TSS = endPos1 

write.csv(df, "~/Wang/output/Paper_figs/Fig5/Peaks_Gpr17_Cor_500KB_Stress_2024-02-12_edit.csv")

