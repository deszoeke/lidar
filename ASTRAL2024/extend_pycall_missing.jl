using PythonCall
const PyObject = Py

# PyObject method interprets Array{Union{T,Missing}} as a
# numpy masked array.
# This allows for plotting with missing values.
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy = pyimport("numpy")
    numpy_ma = numpy.ma
    numpy_ma.array(coalesce.(a,zero(T)), mask=ismissing.(a))
end
