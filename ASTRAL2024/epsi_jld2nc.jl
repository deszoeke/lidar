# save netcdf file with epsilon and se_epsilon for all chunks and heights
using Pkg; Pkg.activate(".")

using NCDatasets
using Dates
# # using Printf
using JLD2

pd = permutedims

# set up data
# input JLD file
data = load("epsilon_data/epsi_stare_vn_chunks_r20260502_2345.jld")
data isa Dict
# add se_epsilon to data dict for saving
data["se_epsi"] = load("epsilon_data/epsi_se.jld")["s_se"]

vars_1d = filter( k -> (data[k] isa Vector), keys(data))
vars_2d = filter( k -> (data[k] isa Matrix), keys(data))
# for k in keys(data)
#     v = data[k]
#     println("$(k) $(ndims(v)) $(size(v)), $(eltype(v))")
# end
# time, height
# for k in vars_1d
#     println("$k $(typeof(data[k])) $(size(data[k]))")
# end
# # time, height
# for k in vars_2d
#     println("$k $(typeof(data[k])) $(size(data[k]))")
# end

nheight = size(data["A"], 2)
nchunk = size(data["A"], 1)


# write the netcdf epsilon file
Dataset("epsilon_data/epsi_stare_vn_chunks_r20260502_2345.nc", "c") do ds
    # dimensions
    defDim(ds, "time", Inf) # record dimension for chunks
    defDim(ds, "height", nheight)

    # height variable
    v_h = defVar(ds, "height", Float64, ("height",))
    v_h[:] = data["height"]
    v_h.attrib["units"] = "meters"

    # time variables: dtime_st and dtime_en
    # Convert DateTime → CF seconds since epoch manually (NCDatasets won't
    # auto-convert when the declared variable type is Float64).
    epoch = DateTime(2024, 4, 1, 0, 0, 0)
    to_cf(v::Vector{DateTime}) = Float64.(Dates.value.(v .- epoch)) ./ 1000.0  # ms → s
    time_attribs = ["units" => "seconds since 2024-04-01 00:00:00 UTC",
                    "calendar" => "proleptic_gregorian"]

    vname = "dtime_st"
    v_t = defVar(ds, vname, Float64, ("time",); attrib=time_attribs)
    v_t.attrib["long_name"] = "start time of chunk"
    v_t[1:nchunk] = to_cf(data[vname])

    vname = "dtime_en"
    v_t = defVar(ds, vname, Float64, ("time",); attrib=time_attribs)
    v_t.attrib["long_name"] = "end time of chunk"
    v_t[1:nchunk] = to_cf(data[vname])

    # time variables
    # for vname in vars_1d
    #     if startswith(vname, "dtime")
    #         println("writing netcdf variable $(vname)")
    #         ds[vname] = defVar(ds, vname, Float64, ("time",)) # CF-compliant Float64 time variable as seconds since epoch
    #         ds[vname].attrib["units"] = "seconds since 2024-04-01 00:00:00 UTC"
    #         ds[vname][:] = data[vname]
    #     end
    # end

    # time-height variables: dimension declaration, data, and attributes
    for vname in vars_2d
        println("writing netcdf variable $(vname)")
        to_f64(x) = Float64.(coalesce.(x, NaN))  # Missing → NaN for netcdf _FillValue
        if vname == "RMSE_tmp" # rename RMSE_tmp -> RMSE_fit
            v = defVar(ds, "RMSE_fit", Float64, ("height", "time"); fillvalue=NaN)
            v.attrib["long_name"] = "root mean square error of structure function fit"
            v[:, 1:nchunk] = pd(to_f64(data[vname])) # transpose to (height, time) order for netcdf convention
        else # all other 2D variables
            v = defVar(ds, vname, Float64, ("height", "time"); fillvalue=NaN)
            v[:, 1:nchunk] = pd(to_f64(data[vname])) # transpose to (height, time) order for netcdf convention
            # attributes
            if contains(vname, "A")
                if startswith(vname, "se_")
                    v.attrib["long_name"] = "standard error of slope A"
                    v.attrib["description"] = "uncertainty of slope A from structure function fit"
                    v.attrib["units"] = "m^(2/3) s^(-2)"
                else
                    v.attrib["long_name"] = "slope of structure function D2 vs. separation^(2/3)"
                    v.attrib["description"] = "A = C_k * epsilon^(2/3)"
                    v.attrib["units"] = "m^(2/3) s^(-2)"
                end
            elseif contains(vname, "epsi")
                if startswith(vname, "se_")
                    v.attrib["long_name"] = "standard error of epsilon"
                    v.attrib["description"] = "uncertainty of TKE dissipation rate from structure function fit"
                    v.attrib["units"] = "m^2 s^(-3)"
                else
                    v.attrib["long_name"] = "TKE dissipation rate"
                    v.attrib["description"] = "epsilon computed from slope A of structure function fit"
                    v.attrib["equation"] = "D2 = Ck * epsilon^(2/3) * separation^(2/3)"
                    v.attrib["units"] = "m^2 s^(-3)"
                end
            elseif vname == "nbins"
                v.attrib["long_name"] = "number of bins in structure function fit"
                v.attrib["units"] = "count"
            end
        end
    end

end

