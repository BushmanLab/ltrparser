for f in analysis_data/GTSP3*1_u[35].R2.fastq.gz; do f2=${f/-1/}; f2=${f2/analysis_data/subsample}; zcat $f | head -n 40000 | gzip > $f2; done
for f in analysis_data/GTSP3*2_u[35].R2.fastq.gz; do f2=${f/-2/}; f2=${f2/analysis_data/subsample}; zcat $f | head -n 40000 | gzip >> $f2; done
for f in analysis_data/GTSP3*3_u[35].R2.fastq.gz; do f2=${f/-3/}; f2=${f2/analysis_data/subsample}; zcat $f | head -n 40000 | gzip >> $f2; done
for f in analysis_data/GTSP3*4_u[35].R2.fastq.gz; do f2=${f/-4/}; f2=${f2/analysis_data/subsample}; zcat $f | head -n 40000 | gzip >> $f2; done

