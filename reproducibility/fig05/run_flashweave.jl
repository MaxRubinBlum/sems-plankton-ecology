using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using FlashWeave
using CSV
using DataFrames

# Run from reproducibility/fig05/. Input tables are expected in data/flashweave/.
repo = normpath(joinpath(@__DIR__, "..", ".."))
inp = joinpath(repo, "data", "flashweave")
out = joinpath(repo, "results", "figure5_cross_domain_associations", "flashweave")
mkpath(out)

feature_path = joinpath(inp, "features_complete_cases.tsv")
meta_path = joinpath(inp, "metadata_complete_cases.tsv")

features = CSV.read(feature_path, DataFrame; delim='\t')
metadata = CSV.read(meta_path, DataFrame; delim='\t')
@assert String.(features[:,1]) == String.(metadata[:,1]) "Feature and metadata sample order differs."

println("Samples: ", nrow(features))
println("Features: ", ncol(features)-1)
println("Metadata variables: ", ncol(metadata)-1)

println("Running FlashWeave-Sensitive with metadata...")
net_meta = learn_network(feature_path, meta_path;
    sensitive=true,
    heterogeneous=false,
    max_k=3,
    alpha=0.01,
    FDR=true,
    verbose=true)
save_network(joinpath(out, "flashweave_sensitive_with_metadata.edgelist"), net_meta)
save_network(joinpath(out, "flashweave_sensitive_with_metadata.gml"), net_meta)

println("Running FlashWeave-Sensitive without metadata...")
net_nometa = learn_network(feature_path;
    sensitive=true,
    heterogeneous=false,
    max_k=3,
    alpha=0.01,
    FDR=true,
    verbose=true)
save_network(joinpath(out, "flashweave_sensitive_without_metadata.edgelist"), net_nometa)
save_network(joinpath(out, "flashweave_sensitive_without_metadata.gml"), net_nometa)

println("FLASHWEAVE RUN COMPLETED SUCCESSFULLY")
