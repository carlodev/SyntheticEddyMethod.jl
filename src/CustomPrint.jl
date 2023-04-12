function Base.show(io::IO, obj::Union{SemEddy,VirtualBox}) 
    get_type = typeof(obj)
    print(io, "$(get_type)\n") 
    for fname in fieldnames(get_type)
        fval = getfield(obj,fname)
        print(io, "$fname = $fval\n")    
    end
end
