cd data
awk '{print $1}' meta_data.csv > sample_id.csv
sed 's/.* //g' meta_data.csv | sed 's/^/,/g' > IBD_val.csv
paste sample_id.csv IBD_val.csv | sed 's/biopsy/biopsy,/g' | sed 's/stool/stool,/g' | sed 's/\s\+//g' > meta_data_filtered.csv
sed -I ’s/#SampleID.*/#SampleID,type_sample,IBD/g' meta_data_filtered.csv
sed -i ''  's/#SampleID.*/#SampleID,type_sample,IBD/g' meta_data_filtered.csv
rm IBD_val.csv sample_id.csv
cd ../ 
