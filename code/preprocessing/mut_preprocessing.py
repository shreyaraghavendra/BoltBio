import os
import gzip
from typing import Counter
import pandas as pd
from utils import data_dir_windows

#root_dir = '../../data/pancancer/TCGA/mutation/download_new/'
#mut_data = data_dir_windows + "mutation"
def combine_maf(combine_file_path = "/combined_maf.maf", data_path = data_dir_windows + "mutation"):
    count = 0
    maf_all = None
    col_names = None
    for subdir, d, files in os.walk(data_path):  
        for fname in files:
            print(fname)
            if fname.endswith('.maf.gz'):
                p = os.path.join(subdir, fname)
                if maf_all is None:
                    maf_all = pd.read_csv(p, sep='\t', comment='#', header=0)
                    col_names = maf_all.columns
                else:
                    print ("[{}] Got {} mutations and {} common columns".format(count, maf_all.shape[0], maf_all.shape[1]))
                    try:
                        maf_new = pd.read_csv(p, sep='\t', comment='#', header=0)
                    except:
                        print ("Exception occured while reading {}... Skipping!".format(p))
                        continue
                    if {'Hugo_Symbol', 'Variant_Classification', 'Tumor_Sample_Barcode'}.issubset(maf_new.columns):
                        maf_join = pd.concat([maf_all, maf_new], ignore_index=True, join='inner')
                        maf_join.dropna(axis=1, inplace=True, how='all') # remove stupid all-na columns
                        if len(maf_join.columns) < len(col_names):
                            print ("Lost {} columns reading {}".format(len(col_names)-maf_join.shape[1], p))
                            print ("Lost columns: {}".format([c for c in col_names if not c in maf_join.columns]))
                            print ("New Maf cols: {}".format(maf_new.columns))
                            print ("Old MAF cols: {}".format(maf_all.columns))
                            print ("Join MAF cols: {}".format(maf_join.columns))
                            col_names = maf_all.columns
                        maf_all = maf_join
                    else:
                        print ("File contains not all columns needed or is empty {}".format(p))
                count += 1

    maf_all.to_csv(combine_file_path, sep='\t',index=False)

def convert_to_annotation_input(input_file, output_file):
    maf_file = pd.read_table(input_file)
    oc_input_df = maf_file[["Chromosome","Start_Position","Strand","Reference_Allele","Tumor_Seq_Allele2"]]
    oc_input_df.to_csv(output_file, header=None, index=None, sep=' ', mode='a')
    

    
combine_maf(combine_file_path = data_dir_windows + "mutation" + '/combined_maf.maf',data_path=data_dir_windows + "mutation")
convert_to_annotation_input(input_file = data_dir_windows + "mutation/combined_maf.maf", output_file= data_dir_windows + "mutation/cravat_input")