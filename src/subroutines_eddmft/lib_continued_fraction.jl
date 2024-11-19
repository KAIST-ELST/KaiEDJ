using LinearAlgebra
using SparseArrays
using ThreadedSparseCSR
using SparseMatricesCSR
export lanczos_algorithm
export continued_fraction

function lanczos_algorithm_orig(H::SparseMatrixCSC{T}, v0::Vector{T}, m::Int) where T
    n = size(H, 1)  # Dimension of the Hamiltonian matrix
    V = Matrix{T}(undef, n, m + 1)  # Matrix to store Lanczos vectors
    T_mat = Matrix{T}(undef, m, m)  # Tridiagonal matrix

    V[:, 1] = normalize(v0)  # Normalize the initial vector

    beta = 0.0
    for j in 1:m
        w       = H * V[:, j]
        alpha   = V[:, j]' * w
        w       -= alpha * V[:, j] .+ (j > 1 ? beta * V[:, j - 1] : zero(T))
        beta    = norm(w)

        T_mat[j, j] = alpha
        if beta < 1e-10
            break
        end
        V[:, j + 1] = w / beta
        # println( "Tjj : ", T_mat[j,j]j)
        if j > 1
            T_mat[j - 1, j] = beta
            T_mat[j, j - 1] = beta
            # println( "Tij : ", T_mat[j,j], T_mat[j-1,j], T_mat[j,j-1] )
        end
    end

    return T_mat, V
end


function lanczos_algorithm(H::SparseMatrixCSC{T}, v0::Vector{T}, m::Int) where T
    n = size(H, 1)  # Dimension of the Hamiltonian matrix
    V = Matrix{T}(undef, n, m + 1)  # Matrix to store Lanczos vectors
    V       = zero(V)
    T_mat   = Matrix{T}(undef, m, m)  # Tridiagonal matrix
    T_mat   = zero(T_mat)

    V[:, 1] = normalize(v0)  # Normalize the initial vector

    # ShowProbArrNonzero( [ v0 ] ) 
    # @show V[:, 1]' * H * V[:,1]

    beta = 0.0
    for j in 1:m
        Hvj       = H * V[:, j]
        alpha   = V[:, j]' * Hvj
        T_mat[j, j] = alpha

        tilde_phi_new   = Hvj - alpha * V[:, j] 
        if j > 1 
            tilde_phi_new  -= beta * V[:, j - 1] 
            T_mat[j-1, j] = beta
            T_mat[j, j-1] = conj(beta)
        end

        beta = norm(tilde_phi_new)
        if beta < 1e-12
            # V[:, j+1] = normalize!(randn(T,n))
            # vnew  = V[:, j+1]
            # for k in 1:j
            #     vnew -= (V[:,k]' * vnew) * V[:,k]
            # end
            # norm_vnew   = norm(vnew)
            # normalize!(vnew)
            # if norm_vnew < 1e-12
            #     break
            # end
            # V[:, j+1]   = vnew
            # beta        = V[:,j]' * tmul(H, vnew)
            break
        else
            V[:, j+1] = normalize(tilde_phi_new)
        end

        # println( "Tjj $j : ", T_mat[j,j])
        # if j > 1
        #     println( format("Tij Tj-1j Tjj-1 : {} {} {} ", T_mat[j,j], T_mat[j-1,j], T_mat[j,j-1]) )
        # end
        # println("")
    end
    # if isnan(sum(T_mat)) 
    #     @show beta
    #     norm_arr    = [ norm(V[:,j]) for j=1:m+1 ]
    #     @show norm_arr
    #     error( "isnan(sum(T_mat))  = true " ) 
    # end

    return T_mat, V
end

function lanczos_algorithm(H::SparseMatrixCSR, v0::Vector{T}, m::Int) where T
    n = size(H, 1)  # Dimension of the Hamiltonian matrix
    V       = Matrix{T}(undef, n, m + 1)  # Matrix to store Lanczos vectors
    V       = zero(V)
    T_mat   = Matrix{T}(undef, m, m)  # Tridiagonal matrix
    T_mat   = zero(T_mat)

    V[:, 1] = normalize(v0)  # Normalize the initial vector

    # ShowProbArrNonzero( [ v0 ] ) 
    # @show V[:, 1]' * H * V[:,1]

    ThreadedSparseCSR.multithread_matmul(BaseThreads())

    beta = 0.0
    for j in 1:m
        vj      = V[:,j]
        Hvj     = tmul(H, vj)   ## H * V[:, j]
        # Hvj      = zero( vj )
        # tmul!( $Hvj, $H, $vj )   ## H * V[:, j]
        alpha_j     = vj' * Hvj
        T_mat[j, j] = alpha_j

        tilde_phi_new   = Hvj - alpha_j * vj
        if j > 1 
            tilde_phi_new   -= beta * V[:, j - 1] 
            T_mat[j-1, j]   = beta
            T_mat[j, j-1]   = conj(beta)
        end

        beta = norm(tilde_phi_new)
        if beta < 1e-12
            # V[:, j+1] = normalize!(randn(T,n))
            # vnew  = V[:, j+1]
            # for k in 1:j
            #     vnew -= (V[:,k]' * vnew) * V[:,k]
            # end
            # norm_vnew   = norm(vnew)
            # normalize!(vnew)
            # if norm_vnew < 1e-12
            #     for ii in j+1:m
            #         V[ii, ii] = 1.0
            #     end
            #     break
            # end
            # V[:, j+1]   = vnew
            # beta        = V[:,j]' * tmul(H, vnew)
            break
        else
            V[:, j+1] = norm(tilde_phi_new)<1e-12 ? zero(tilde_phi_new) :  normalize(tilde_phi_new)
        end

        # println( "Tjj $j : ", T_mat[j,j])
        # if j > 1
        #     println( format("Tij Tj-1j Tjj-1 : {} {} {} ", T_mat[j,j], T_mat[j-1,j], T_mat[j,j-1]) )
        # end
        # println("")
    end
    # if isnan(sum(T_mat)) 
    #     @show beta
    #     norm_arr    = [ norm(V[:,j]) for j=1:m+1 ]
    #     @show norm_arr
    #     error( "isnan(sum(T_mat))  = true " ) 
    # end

    return T_mat, V
end

# function continued_fraction(M::Matrix{T}, m::Int, energy::T) where T
#     G = zero(T)
#     M[1, 1] += energy
#     M[1, 2] += 1 / M[1, 1]
#     for j in 2:m-1
#         M[j, j] += energy - M[j - 1, j] * M[j, j - 1]
#         M[j, j + 1] += 1 / M[j, j]
#         G = M[j, j + 1]
#     end
#     return G
# end

function continued_fraction(TriMat::Matrix{T}, m::Int, energy::T) where T
    M   = deepcopy(TriMat)
    G   = zero(T)
    # println( "Mmm : ", M[m,m], M[m,m-1] )
    M[m, m  ] += energy
    M[m, m-1] *= 1 / M[m, m]
    # println( "Mmm : ", M[m,m], M[m,m-1] )
    # println("(m-1):-1:2")
    # println( (m-1):-1:2 )
    # println( collect((m-1):-1:2) )
    for j in (m-1):-1:2
        # println( "Mjj : ", M[j,j], M[j+1,j], M[j,j+1] )
        M[j, j      ] += energy - M[j + 1, j] * M[j, j + 1]
        M[j, j - 1  ] *= 1 / M[j, j]
    end
    M[1, 1] += energy - M[2,1] * M[1,2]
    G   = 1. / M[1,1]
    return G
end

# # Example usage
# # Define the Hamiltonian matrix (sparse matrix) and the initial vector
# # H = ... (define your Hamiltonian matrix)
# # v0 = ... (define your initial vector)
# 
# # Set the number of Lanczos iterations (m)
# m = 10
# 
# # Run the Lanczos algorithm
# T, V = lanczos_algorithm(H, v0, m)
# @show typeof(T)
# @show size(T)
# @show typeof(V)
# @show size(V)
# 
# # Choose an energy value for the Green's function
# energy = 0.0  # Replace this with your desired energy value
# 
# # Calculate the continued fraction representation of the Green's function
# G = continued_fraction(T, m, energy)
# 
# println("Green's function at energy ", energy, ": ", G)
# 
