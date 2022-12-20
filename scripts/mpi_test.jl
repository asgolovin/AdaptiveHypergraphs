using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
println("Hello world, I am $(rank) of $(size)")
MPI.Barrier(comm)

secret = nothing
if rank == 0
    secret = 42
end

secret = MPI.bcast(secret, 0, comm)

println("The secret is $secret from rank $rank")

MPI.Finalize()