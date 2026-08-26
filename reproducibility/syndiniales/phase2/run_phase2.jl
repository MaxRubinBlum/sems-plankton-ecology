using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using FlashWeave

length(ARGS) == 2 || error("Usage: julia --project=. run_phase2.jl INPUT_DIR OUTPUT_DIR")
inp = abspath(ARGS[1]); out = abspath(ARGS[2]); mkpath(out)
for reg in ["upper_le_220m", "deep_gt_220m"]
    f = joinpath(inp, "features_" * reg * ".tsv")
    m = joinpath(inp, "metadata_" * reg * ".tsv")
    println("Running FlashWeave-Sensitive: ", reg)
    net = learn_network(f, m; sensitive=true, heterogeneous=false, max_k=3, alpha=0.01, FDR=true, verbose=true)
    save_network(joinpath(out, "phase2_" * reg * ".edgelist"), net)
    save_network(joinpath(out, "phase2_" * reg * ".gml"), net)
end
