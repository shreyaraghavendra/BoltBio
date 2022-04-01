import pandas as pd

'''Normalize all the annotation to the range 0 to 1 where 0 is benign and 1 is damaging'''
def normalize_annotation(input_file):
    ann_file = pd.read_excel(input_file,sheet_name="Variant",header=[0,1])
    print(ann_file["SIFT"])
    print(ann_file.head())



normalize_annotation("data\\mutation\\annotated_mutation.xlsx")
