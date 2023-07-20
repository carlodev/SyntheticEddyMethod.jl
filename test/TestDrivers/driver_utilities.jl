
x = rand(5)
y = rand(7)
z = rand(17)
#Computing the velocity in the middle of the VirtualBox domain
vector_points = create_vector_points(x, y, z)

@test length(vector_points) == length(x) * length(y) * length(z)



@test isapprox(compute_Ek([1.2, 0.1, 0.3], 1.0), (0.2^2 + 0.1^2 + 0.3^2) * 0.5)
