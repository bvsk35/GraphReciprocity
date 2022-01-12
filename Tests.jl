
#cd("C:\\Users\\dcirk\\Documents\\CSCE689\\Project")

##### Facebook Analysis ######

# Load in Facebook Data
using DelimitedFiles, MatrixMarket, JSON
include("MaxRecipFunctions.jl")
edgelist = readdlm("Datasets\\facebook-wall.txt")
edgelist = convert(Matrix{Int64}, edgelist[:, 1:2])

# Convert to adjacency matrix
n_edges = size(edgelist)[1]
n_nodes = maximum(edgelist)
A = sparse(edgelist[:, 1], edgelist[:, 2], repeat([1], n_edges), n_nodes, n_nodes)
# Rewire algorithm
time_facebook = @timed begin
A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
facebook_recip = recip/n_edges
# reciprocity of rewired graph
facebook_rewire_recip = recip_rewire/n_edges
# upper bound
facebook_upper_bound = upperbound/n_edges
facebook_edges = n_edges
facebook_nodes = n_nodes

##### EU Email Analysis ######

A = MatrixMarket.mmread("Datasets\\email-Eu-core-temporal\\email-Eu-core-temporal.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]
# Rewire algorithm
time_email_EU = @timed begin
A_rewire = greedyRewire(A)
end

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
email_EU_recip = recip/n_edges
# reciprocity of rewired graph
email_EU_rewire_recip = recip_rewire/n_edges
# upper bound
email_EU_upper_bound = upperbound/n_edges
email_EU_edges = n_edges
email_EU_nodes = n_nodes

##### Slashdot Analysis ######

A = MatrixMarket.mmread("Datasets\\soc-Slashdot0811\\soc-Slashdot0811.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]
# Rewire algorithm
time_slashdot = @timed begin
    A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
slashdot_recip = recip/n_edges
# reciprocity of rewired graph
slashdot_rewire_recip = recip_rewire/n_edges
# upper bound
slashdot_upper_bound = upperbound/n_edges
slashdot_edges = n_edges
slashdot_nodes = n_nodes

#### Mathoverflow Analysis ####

A = MatrixMarket.mmread("Datasets\\sx-mathoverflow\\sx-mathoverflow.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]
# Rewire algorithm
time_mathoverflow = @timed begin
    A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
mathoverflow_recip = recip/n_edges
# reciprocity of rewired graph
mathoverflow_rewire_recip = recip_rewire/n_edges
# upper bound
mathoverflow_upper_bound = upperbound/n_edges
mathoverflow_edges = n_edges
mathoverflow_nodes = n_nodes

#### College Message Analysis ####

A = MatrixMarket.mmread("Datasets\\CollegeMsg\\CollegeMsg.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]
# Rewire algorithm
time_collegemsg = @timed begin
    A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
collegemsg_recip = recip/n_edges
# reciprocity of rewired graph
collegemsg_rewire_recip = recip_rewire/n_edges
# upper bound
collegemsg_upper_bound = upperbound/n_edges
collegemsg_edges = n_edges
collegemsg_nodes = n_nodes

#### Epinions Analysis ####
A = MatrixMarket.mmread("Datasets\\soc-Epinions1\\soc-Epinions1.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]

# Rewire algorithm
time_epinions = @timed begin
    A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
epinions_recip = recip/n_edges
# reciprocity of rewired graph
epinions_rewire_recip = recip_rewire/n_edges
# upper bound
epinions_upper_bound = upperbound/n_edges
epinions_edges = n_edges
epinions_nodes = n_nodes

#### Twitter Analysis ####
A = MatrixMarket.mmread("Datasets\\higgs-twitter\\higgs-twitter_mention.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]

# Rewire algorithm
time_twitter = @timed begin
    A_rewire = greedyRewire(A)
end

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
twitter_recip = recip/n_edges
# reciprocity of rewired graph
twitter_rewire_recip = recip_rewire/n_edges
# upper bound
twitter_upper_bound = upperbound/n_edges
twitter_edges = n_edges
twitter_nodes = n_nodes

#### PA Analysis ####
A = MatrixMarket.mmread("Datasets\\PAdata.mtx")
# Convert to adjacency matrix
n_edges = sum(A)
n_nodes = size(A)[1]

# Rewire algorithm
time_PA = @timed begin
    A_rewire = greedyRewire(A)
end
    

# Reciprocity before and after
recip = compute_reciprocity(A)
recip_rewire = compute_reciprocity(A_rewire)

out_deg = transpose(sum(A, dims = 1))
in_deg = sum(A, dims = 2)
min_deg = zeros(n_nodes)
for i in 1:n_nodes
    min_deg[i] = min(in_deg[i], out_deg[i])
end
upperbound = sum(min_deg)

# reciproicty of original graph
PA_recip = recip/n_edges
# reciprocity of rewired graph
PA_rewire_recip = recip_rewire/n_edges
# upper bound
PA_upper_bound = upperbound/n_edges

PA_edges = n_edges
PA_nodes = n_nodes


output_data =
["Facebook" facebook_recip facebook_rewire_recip facebook_upper_bound facebook_nodes facebook_edges time_facebook[2]
"EU Email" email_EU_recip email_EU_rewire_recip email_EU_upper_bound email_EU_nodes email_EU_edges time_email_EU[2]
"Slashdot" slashdot_recip slashdot_rewire_recip slashdot_upper_bound slashdot_nodes slashdot_edges time_slashdot[2]
"Mathoverflow" mathoverflow_recip mathoverflow_rewire_recip mathoverflow_upper_bound mathoverflow_nodes mathoverflow_edges time_mathoverflow[2]
"College Msg" collegemsg_recip collegemsg_rewire_recip collegemsg_upper_bound collegemsg_nodes collegemsg_edges time_collegemsg[2]
"Epinions" epinions_recip epinions_rewire_recip epinions_upper_bound epinions_nodes epinions_edges time_epinions[2]
"Twitter" twitter_recip twitter_rewire_recip twitter_upper_bound twitter_nodes twitter_edges time_twitter[2]
"Preferential" PA_recip PA_rewire_recip PA_upper_bound PA_nodes PA_edges time_PA[2]] 

#writedlm("rewiringResults.txt", output_data)