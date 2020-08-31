using Random
using LinearAlgebra
using Distributions
using JLD
using DataFrames
using CSV
using DelimitedFiles
using TypedTables
using Query
using MATLAB
using DataFramesMeta
####################### Ques 1 ###########################################
# Part (a)

function q1()
    Random.seed!(1234)
    A = rand(Uniform(-1,2), 10,7)
    B = rand(Normal(-2,15),10,7)
    C = hcat(A[1:5,1:5],B[1:5,6:7])
    D = similar(A)
    for i in 1:10, j = 1:7
        if A[i,j] <= 0
            D[i,j] = A[i,j]
        else
            D[i,j] = 0
        end
    end

    # Part(b)
    println(length(A)) ## To list the number of elements of the matrix
    size(A)   ## To describe the size of the matrix

    # Part (c) : List the number of unique elements in the matrix D #
    length(unique(D))

    # Part(d) : Converting a matrix to a vector
    E = reshape(B,70)  # Using Reshape Function
    E  = B[:]          # Simple Way

    # Part(e): Creating a three-dimensional vector using A in the first column of the 3rd dimension and B in the second column of the 3rd dimension
    F = cat(A,B,dims=3) # Using Concatenation

    # Part(f): Using permutedims to transform F from (10*7*2)  to (2*10*7)
    #F = permutedims(F,(2,10,7)) # Not Working
    F = reshape(F,(2,10,7))     # Reshape Working

    # Part(g): Kronecker Product of B and C 
    G = kron(A,B)
    #kron(C,F)   ## Kron does not work well for 3-D matrices

    ### Part(h): Saving Matrix A,B,C,D,E,F,G as a .jld file called matrixpractice

    # Saving all the matrices with their respective unique identifiers
    save("Solutions_Problem_Sets/matrixpractice.jld",
    "A",A,"B",B,"C",C,"D",D,"E",E,"F",F,"G",G)

    # Fetching Only the Matrix A 
    #load("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/matrixpractice.jld","A")

    ### Part(i): Save only the matrices A, B, C, and D as a .jld file called firstmatrix.
    save("Solutions_Problem_Sets/firstmatrix.jld",
    "A",A,"B",B,"C",C,"D",D)

    ### Part(j): Export C as a .csv file called Cmatrix. You will first need to transform C into a Dataframe
    C = convert(DataFrame,C) # Converting to dataframe

    # Writitng a csv file
    CSV.write("Solutions_Problem_Sets/Cmatrix.csv",C)

    ### Part(h):Export D as a tab-delimited .dat file called Dmatrix. You will first need to transform D into a DataFrame.
    D = convert(DataFrame,D)

    # Writing a tab delimited .dat file
    CSV.write("Solutions_Problem_Sets/Dmatrix.dat", D)
    #CSV.File("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/Dmatrix.dat") |> DataFrame!


    ## Part(k):  Wrapping it in function
    return A,B,C,D
end

A,B,C,D = q1()

############### Ques2: Practice with Loops and Comprehensions ##########################

#=  Part(a):Write a loop or use a comprehension that computes the element-by-element product of
A and B. Name the new matrix AB. Create a matrix called AB2 that accomplishes this
task without a loop or comprehension. =#
function q2(A,B,C,D)
    AB = zeros(10,7)  # Using Loop
    for i in 1:length(A)
        AB[i] = A[i]*B[i]
    end

    AB1 = [ row[1]*row[2] for row in zip(A,B)] # Using Comprehensions

    AB2 = A .* B    # With dot this is an element by element operation, and without dot this will be a dot product.


    #= Part(b): Write a loop that creates a column vector called Cprime which contains only the elements
    of C that are between -5 and 5 (inclusive). Create a vector called Cprime2 which
    does this calculation without a loop. =#

    Cprime = Vector{Float64}()         # First doing this calculation without a loop
    for item in eachrow(C), i in item
        if   -5 <= i <= 5
            append!(Cprime,i)
        end
    end
    Cprime
    ###### DOUBT ###########
    # I am not able to perform this without loop
    #X[ :,2,1]

    #= Part 2(c): Please Refer the Problem Set 1 =#
    X = zeros((15169,6,5))      
    # Using for loop
    for i in 1:size(X)[3], j in 1:size(X)[2]
        if j == 1
            X[:,j,i] .= 1
        elseif  j==2 
            a2 = rand(Binomial(1,(0.75*(6-i))/5),15169) 
            X[:,j,i] = a2
        elseif j==3 
            a3 = rand(Normal(15+i-1,5*(i-1)),15169)
            X[:,j,i] = a3
        elseif j == 4 
            a4 = rand(Normal((pi*(6-i))/3,1/ℯ),15169)
            X[:,j,i] = a4
        elseif  j == 5 
            a5 = rand(Normal(12,2.19),15169) ## Limitation: Not Know how to sample from Discreet Normal R.V.
            X[:,j,i] = a5
        elseif j == 6
            a6 = rand(Binomial(20,0.5),15169) 
            X[:,j,i]  = a6   
        end
    end

    #= Part 2(d):Use comprehensions to create a matrix b which is K T and whose elements evolve
    across time in the following fashion:=#

    #= a =[2;3;4;5;6] Using comprehensions with conditionals
    b =[0;7;8;9;10]
    c = [x < y ? x : y for (x, y) in zip(a, b)] =#

    K = 6
    T = 5
    β = zeros(K,T)
    β[1,:] = [1:0.25:2;]
    β[2,:] = [ log(x) for x in 1:T ]
    β[3,:] = [-sqrt(x) for x in 1:T ]
    β[4,:] = [ exp(x)-exp(x+1) for x in 1:T ]
    β[5,:] = [x for x in 1:T]
    β[6,:] = [x/3 for x in 1:T]

    #= Part 2(e): Use comprehensions to create a matrix Y which is N*T.... =#
    Y = zeros((15169,5,5))

    for i in 1:5
        ϵₜ = rand(Normal(0,0.36),(15169,5))
        Y[:,:,i] = X[:,:,i]*β + ϵₜ
    end

end
q2(A,B,C,D)