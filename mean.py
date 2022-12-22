import pandas

DepthA=pandas.read_csv("Haslea_O/16S/Matam/Depth/.depth",sep='\t')
DepthA.columns = ['contig', 'locus', 'depth']
DepthA=DepthA.drop(columns=['locus'])
MeanA=DepthA.groupby(['contig']).mean()
MeanA.to_csv("Haslea_O/16S/Matam/Depth/.meandepth",sep='\t',header=False)
