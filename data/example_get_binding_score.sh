Rscript Calculate_binding_scores.R Cmai_outputFile_with_rank.csv BCR_metrics.csv 0.03 "(10,15,20)" BCR_binding_score.csv
# Cmai_outputFile_with_rank.csv is the file name of your Cmai results and the column 'Rank' is required.
# BCR_metrics.csv is the file name containing your BCR information with the column BCR_id and bcr_metric which is usually the clone size for each BCR.
# 0.03 is the threshold
# "(10,15,20)" are the cutoffs which the length varies.
# BCR_binding_score.csv is the output file of your binding scores.
