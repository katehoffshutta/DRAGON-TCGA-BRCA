from netZooPy.dragon import *
import pandas
import os
import sys

#expr_file = "../../data/interim/gene_expression_processed.tsv"
#meth_file = "../../data/interim/methylation_processed.tsv"

expr_file = sys.argv[1]
meth_file = sys.argv[2]
out_dir = sys.argv[3]

expr = pandas.read_table(expr_file,header=0,index_col=0)
meth = pandas.read_table(meth_file,header=0,index_col=0)

# which ids are in both?
all_data = pandas.merge(expr,meth,on="TCGA_barcode",how="inner")

# subset only the ones we want for dragon

meth_data = all_data.filter(regex='methylation')
exp_data = all_data.filter(regex='expr')

# run dragon

lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(meth_data,exp_data)
print(lambdas)

r = get_partial_correlation_dragon(meth_data,exp_data,lambdas)
newnames = sum([meth_data.columns.tolist(),exp_data.columns.tolist()],[])

out_dir_long = "data/processed/" + out_dir
if not os.path.exists(out_dir_long):
    os.mkdir(out_dir_long)

df = pandas.DataFrame(r,columns=newnames,index=newnames)
df.to_csv(out_dir_long + "/dragon_mat.tsv",sep="\t")

all_data.to_csv(out_dir_long + "/dragon_input_mat.tsv",sep="\t")

n = exp_data.shape[0]
p1 = exp_data.shape[1]
p2 = meth_data.shape[1]

print(n)
print(p1)
print(p2)

adj_p_vals, p_vals = estimate_p_values_dragon(r, n, p1, p2, lambdas)

df_adj = pandas.DataFrame(adj_p_vals,columns=newnames,index=newnames)
df_adj.to_csv(out_dir_long + "/dragon_adj_p.tsv",sep="\t")

df_raw = pandas.DataFrame(p_vals,columns=newnames,index=newnames)
df_raw.to_csv(out_dir_long + "/dragon_raw_p.tsv",sep="\t")

