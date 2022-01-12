using SparseArrays, LinearAlgebra

# Given a triangle, rewire it. Returns A_leftover and A_reciprocal 
function rewire_triangle(nodes, A_leftover, A_reciprocal)
    # Find smallest number of edges in a triangle side
    min_edges = min(A_leftover[nodes[1], nodes[2]], A_leftover[nodes[2], nodes[3]], A_leftover[nodes[3], nodes[1]])
    n_recip = 0
    # If smallest number of edges 2 or more, do even rewiring
    if min_edges > 1
        # Number of reciprocal edges to add
        n_recip = (min_edges - mod(min_edges, 2))/2
        A_reciprocal[nodes[2], nodes[1]] = A_reciprocal[nodes[1], nodes[2]] = A_reciprocal[nodes[2], nodes[1]] + n_recip
        A_reciprocal[nodes[3], nodes[2]] = A_reciprocal[nodes[2], nodes[3]] = A_reciprocal[nodes[3], nodes[2]] + n_recip
        A_reciprocal[nodes[1], nodes[3]] = A_reciprocal[nodes[3], nodes[1]] = A_reciprocal[nodes[1], nodes[3]] + n_recip
        A_leftover[nodes[1], nodes[2]] = A_leftover[nodes[1], nodes[2]] - 2 * n_recip
        A_leftover[nodes[2], nodes[3]] = A_leftover[nodes[2], nodes[3]] - 2 * n_recip
        A_leftover[nodes[3], nodes[1]] = A_leftover[nodes[3], nodes[1]] - 2 * n_recip
    end
    # If applicable, do the one rewiring
    if (mod(min_edges, 2) != 0)
        A_leftover[nodes[1], nodes[2]] = A_leftover[nodes[1], nodes[2]] - 1
        A_leftover[nodes[2], nodes[3]] = A_leftover[nodes[2], nodes[3]] - 1
        A_leftover[nodes[3], nodes[1]] = A_leftover[nodes[3], nodes[1]] - 1
        A_leftover[nodes[1], nodes[1]] = A_leftover[nodes[1], nodes[1]] + 1
        A_reciprocal[nodes[3], nodes[2]] = A_reciprocal[nodes[2], nodes[3]] = A_reciprocal[nodes[3], nodes[2]] + 1
    end
    
    return A_reciprocal, A_leftover
end


# Function that removes all reciprocal edges from a graph
function clean_reciprocal(A)
    A_reciprocal = spzeros(size(A)[1], size(A)[1])
    for j in 1:size(A)[1]
        for r in nzrange(A, j)
            i = rowvals(A)[r]
            if i > j
                A_reciprocal[i, j] = min(A[i, j], A[j, i])
            end
        end
    end
    A_reciprocal = A_reciprocal + transpose(A_reciprocal)
    A_leftover = A - A_reciprocal

    return A_reciprocal, A_leftover
end

function greedyRewire(A)
    # Find number of nodes and edges
    n_edges = size(A)[1]
    n_nodes = size(A)[1]
    # Clean recoprical edges from graph
    A_reciprocal, A_leftover = clean_reciprocal(A)
    # Loop over all triangles
    for j in 1:size(A)[1]
        # Loop over all nonzero elements in column j
        for r in nzrange(A, j)
            i = rowvals(A)[r]
            if i != j # Don't count self loops
                # Loop over intermediate nodes (i.e. see if there exists a k that connects (j, k) and (k, i))
                for q in nzrange(A, i)
                    k = rowvals(A)[q]
                    if (j != k) & (i != k) & (A_leftover[j, k] > 0) & (A_leftover[k, i] > 0) & (A_leftover[i, j] > 0)
                        # Rewire first found triangle
                        A_reciprocal, A_leftover = rewire_triangle([k, i, j], A_leftover, A_reciprocal)
                    end
                end
            end
        end
        # Print which node you are on
        print(j)
        print("\n")
    end

    # Order self loops by max degree
    A_leftover_diag = A_leftover[diagind(A_leftover)]
    A_diag_sort = sortperm(-A_leftover_diag)
    k = 0
    deg_min = 0
    # Loop over all selfloops
    for i in 1:(size(A_leftover_diag)[1] - 1)
        k = 1
        while A_leftover_diag[A_diag_sort[i]] > 0
            deg_min = min(A_leftover_diag[A_diag_sort[i + k]], A_leftover_diag[A_diag_sort[i]])
            if deg_min == 0 
                break
            end
            A_reciprocal[A_diag_sort[i], A_diag_sort[i + k]] = A_reciprocal[A_diag_sort[i], A_diag_sort[i + k]] + deg_min
            A_reciprocal[A_diag_sort[i + k], A_diag_sort[i]] = A_reciprocal[A_diag_sort[i + k], A_diag_sort[i]] + deg_min
            A_leftover[A_diag_sort[i], A_diag_sort[i]] = A_leftover[A_diag_sort[i], A_diag_sort[i]] - deg_min
            A_leftover[A_diag_sort[i + k], A_diag_sort[i + k]] = A_leftover[A_diag_sort[i + k], A_diag_sort[i + k]] - deg_min
            A_leftover_diag[A_diag_sort[i + k]] = A_leftover_diag[A_diag_sort[i + k]] - deg_min 
            A_leftover_diag[A_diag_sort[i]] = A_leftover_diag[A_diag_sort[i]] - deg_min
            k = k + 1
        end
        print(i)
        print("\n")
    end

    A_rewire = A_reciprocal + A_leftover
    return A_rewire
end


# Compute number of reciprocal edges
function compute_reciprocity(A)
    A_copy = spzeros(size(A)[1], size(A)[1])
    i = 0
    for j in 1:size(A, 2)
        for r in nzrange(A, j)
            i = rowvals(A)[r]
            A_copy[i, j] = min(A[i, j], A[j, i])
        end
    end
    reciprocal_edges = 2 * sum(triu(A_copy, 1))
    return reciprocal_edges
end

