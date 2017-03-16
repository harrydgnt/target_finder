#HEATMAP 


work_dir <- "/Users/harryyang/Documents/Research/Class/Com\ Sci\ 225/target_finder/"
setwd(work_dir)

input_matrix <- paste0(work_dir, "T1D_binarized_H3K4Me1_pos_state_matrix.txt")
mat <- data.matrix(read.table(input_matrix, header = TRUE, sep = "\t", row.names = 1 ))
mat <- mat[order(rownames(mat)),]

require(ggplot2)
require(pheatmap)
require(RColorBrewer)
RColorBrewer::display.brewer.all()

getPalette = colorRampPalette(brewer.pal(9,"RdYlBu"))

# T1D_heatmap <- heatmap(mat, Rowv = NA, col = cm.colors(15), scale="column")
pheatmap(t(mat), color = rev(getPalette(15)), width=40, height = 10, filename = paste0(work_dir, "T1D_Heatmap_binary_H3K4me1.png"))

