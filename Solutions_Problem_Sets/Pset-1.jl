using Random
using LinearAlgebra
using Distributions
using JLD
using DataFrames
using CSV
using DelimitedFiles
using TypedTables
####################### Ques 1 ###########################################
# Part (a)
a=1
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
    save("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/matrixpractice.jld",
    "A",A,"B",B,"C",C,"D",D,"E",E,"F",F,"G",G)

    # Fetching Only the Matrix A 
    #load("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/matrixpractice.jld","A")

    ### Part(i): Save only the matrices A, B, C, and D as a .jld file called firstmatrix.
    save("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/firstmatrix.jld",
    "A",A,"B",B,"C",C,"D",D)

    ### Part(j): Export C as a .csv file called Cmatrix. You will first need to transform C into a Dataframe
    C = convert(DataFrame,C) # Converting to dataframe

    # Writitng a csv file
    CSV.write("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/Cmatrix.csv",C)

    ### Part(h):Export D as a tab-delimited .dat file called Dmatrix. You will first need to transform D into a DataFrame.
    D = convert(DataFrame,D)

    # Writing a tab delimited .dat file
    CSV.write("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/Dmatrix.dat", D)
    #CSV.File("F:/R_Project_Directories/Trial_PROJECT/Structural_Modelling_2020/Solutions_Problem_Sets/Dmatrix.dat") |> DataFrame!


    ## Part(k):  Wrapping it in function
    return A,B,C,D
end

A,B,C,D = q1()
