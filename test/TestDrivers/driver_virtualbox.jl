
σ = 0.1
b = 5.0
a = 0.0

#Defining the Virtual Box domain
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y, z, σ)

@test Vboxinfo.σ == [σ, σ, σ]
@test Vboxinfo.X_end == σ
@test Vboxinfo.X_start == -σ
@test Vboxinfo.Y_end == b + σ
@test Vboxinfo.Y_start == a - σ
@test Vboxinfo.Z_end == b + σ
@test Vboxinfo.Z_start == a - σ

σv = [0.1, 0.05, 0.02]
σx, σy, σz = σv


Vboxinfo = VirtualBox([0.0],y, z, σv)

@test Vboxinfo.σ == [σx, σy, σz]
@test Vboxinfo.X_end == σx
@test Vboxinfo.X_start == -σx
@test Vboxinfo.Y_end == b + σy
@test Vboxinfo.Y_start == a - σy
@test Vboxinfo.Z_end == b + σz
@test Vboxinfo.Z_start == a - σz


