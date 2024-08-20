# module timing_lidar: methods for syncing lidar, VectoNav, POSMV clocks and motions

module timing_lidar

##for methods
using Dates
using Statistics
using StatsBase
using Interpolations
##methods not used
# using SignalAlignment
# using Optim

##data IO and visualization
# using JLD2
# using MAT
# using PyPlot

export grid_search

# helper functions
pd = permutedims
m2n(x) = ismissing(x) ? NaN : x
anom(x; dims=1) = x .- mean(x; dims=dims)


## vectornav timing functions

# Define the GPS epoch
const GPS_EPOCH = DateTime(1980, 1, 6) # DateTime
# const GPS_OFFSET = Dates.datetime2epochms(GPS_EPOCH) # Integer time interval since DateTime epoch

# Function to calculate the number of leap seconds between two dates
function leap_seconds(date::DateTime)
    leap_seconds_list = [
        DateTime(1981, 7, 1), DateTime(1982, 7, 1), DateTime(1983, 7, 1),
        DateTime(1985, 7, 1), DateTime(1988, 1, 1), DateTime(1990, 1, 1),
        DateTime(1991, 1, 1), DateTime(1992, 7, 1), DateTime(1993, 7, 1),
        DateTime(1994, 7, 1), DateTime(1996, 1, 1), DateTime(1997, 7, 1),
        DateTime(1999, 1, 1), DateTime(2006, 1, 1), DateTime(2009, 1, 1),
        DateTime(2012, 7, 1), DateTime(2015, 7, 1), DateTime(2017, 1, 1)
    ]
    return count(ls -> ls <= date, leap_seconds_list)
end

# "convert gpstime (nanoseconds) to DateTime"
# function gps2dt(gpstime::Integer)
#     t_d = floor(gpstime / (86_400 * 1_000_000_000)) # nanosecond -> integer day
#     t_ms = round(Int64, (gpstime % (86_400 * 1_000_000_000) ) / 1_000_000 - 1_000*leap_seconds(t_d)) # nansecond -> millisecond of day
#     GPS_EPOCH + Day(t_d) + Millisecond(t_ms)
# end

"convert gpstime (nanoseconds) to DateTime."
function gps2dt(gpstime::Integer)
    # leap seconds not needed when using time elapsed
    # Millisecond(Integer) -> millisecond timedelta
    GPS_EPOCH + Millisecond(round(Int64, gpstime / 1_000_000))
end

"convert gpstime (nanoseconds) to DateTime epoch milliseconds (integer)."
function gps2ms(gpstime::Integer)
    round(Int64, gpstime / 1_000_000 ) + Dates.datetime2epochms(GPS_EPOCH)
end

dt2gpsns(dt) = Dates.value(Millisecond( dt-GPS_EPOCH )) * 1_000_000 # -> Integer
GPS_MORETHAN = dt2gpsns(DateTime(2024,1,1))

# tests
Dates.datetime2epochms(DateTime(0)) # 0 OK
Dates.epochms2datetime( Dates.datetime2epochms(GPS_EPOCH) ) # 1980-01-06T00:00:00 OK
gps2dt( dt2gpsns(DateTime(2024,4,29,5,25,42,128)) ) # 2024-04-29T05:25:42.128 OK
DateTime(0) + Millisecond(gps2ms( dt2gpsns(DateTime(2024,4,29,5,25,42,128)) )) # 2024-04-29T05:25:42.128 OK
DateTime(0) + Millisecond( 63881587542128 ) # # 2024-04-29T05:25:42.128 OK
gps2ms( dt2gpsns(DateTime(2024,4,29,5,25,42,128)) ) # 63881587542128 OK

"""
itp(dtq) = interpolatedt(dtk, xk, Gridded(Linear()))(Dates.datetime2epochms(dtq))
interpolate using datetimes as coordinates
"""
interpolatedt(dt, x...) = interpolate(Dates.datetime2epochms(dt), x...)

"short wrapper for Dates.value(Millisecond(t))"
DVM(t::TimePeriod) = Dates.value(Millisecond(t))
DVM(t::Integer) = Dates.value(Millisecond(t))
DVM(t::Real) = DVM(round(Int64,t))
DVM(t::AbstractVector{T}) where T<:TimePeriod = DVM.(t)
function DVM(dt::DateTime)
    error("do not use DVM(::DateTime)")
end

"""
Millisecond( round(Integer, v) )
Convert a number to a Millisecond, rounding if necessary. Quasi-inverse of DVM.
"""
MRi(v::Number) = Millisecond( round(Integer, v) ) # convert a number to a Millisecond, rounding if necessary

"""
Interpolate between imprecisely recorded jumps in time.
Need Float64 precision for interpolation
"""
function intbetweentime(prectime::Vector{Float64})
    # monotonic, OK
    ii0 = map(x-> x<=0, diff(prectime)) 
    # ind0 = boolean 1 for imprecise times that should be interpolated
    ind0 = [0; findall(ii0)] .+ 1   # interpolate to these, increment because of diff()
    ind1 = findall(.!ii0) .+ 1 # interpolate from these

    # interpolate (and extrapolate) the inbetween times
    itp = extrapolate( 
            interpolate((ind1,), prectime[ind1], Gridded(Linear())), 
            Line() )
    # fltind0 = filter(i-> i>ind1[1] && i<ind1[end], ind0)
    # prectime[fltind0] = itp(fltind0)
    prectime[ind0] .= itp(ind0)
    # or just just reinterp all indices
    # prectime = itp(axes(prectime))
    return prectime
end


# methods for averaging DateTimes and TimePeriod

"meandt(dt, ...)"
function meandt(dt::AbstractVector{DateTime}, args...)
    ms = @. Dates.value(Millisecond(dt - dt[1]))
    return dt[1] .+ Millisecond(round(Int64, mean(ms, args...)))
end
# "mean for datetimes" updated 2024-08-01
# meandt(dt::AbstractVector{DateTime}) = first(dt) + MRi( mean(DVM.(dt .- first(dt))) )
"meandt(t, ...)"
function meandt(dt::AbstractVector{T}, args...) where T<:TimePeriod
    ms = @. Dates.value(Millisecond(dt))
    return Millisecond(round(Int64, mean(ms, args...)))
end

# find regressing GPS times (those that go backward)
iregress(gpstime) = findall(x-> x<0, diff(gpstime))
# ba(ii) = map( i->i.+[0,1], ii)[:] # ind before and after backward step
ba(ii) = (ii.+[0; 1])[:]
dtregress(gpstime) = gps2dt.( gpstime[ ba(iregress(gpstime)) ] )

"interpolate and add offsets to DateTime of precise time from vectornav."
function precise_dt(vn_raw_dt, gpstime)
    RASP_BASE_TIME = floor(vn_raw_dt[1], Dates.Day(1)) # totally bogus unset offset clock time
    # need Float64 precision for interpolation
    # RASP_BASE_TIME is a datetime before the start of the data.
    
    # add offsets
    prectime = intbetweentime( Float64.(DVM( vn_raw_dt .- RASP_BASE_TIME )) )
    precdt_vn_clock = RASP_BASE_TIME .+ Dates.Millisecond.(round.(Int64, prectime))
    offset_gps_vn = gps2dt(gpstime[1]) - precdt_vn_clock[1] # OK for leg 1
    precdt = precdt_vn_clock + offset_gps_vn # move to absolute GPS clock
    return precdt
    # next linearly stretch to minimize bias from GPS
end

"""
Stretch VectorNav millisecond-interpolated
time linearly to match best GPS times at 
the start and 0-1 minute after the reset.
"""
function stretch_vn_to_gps(precdt, gpstime)
    i0 = 2 # step starting good GPS times
    ir = iregress(gpstime)[1]+1 # step after GPS reset
    # delta = meandt(gps2dt.(gpstime[ir:(ir+1200)])) - meandt(precdt[ir:(ir+1200)])
    delta = meandt(gps2dt.(gpstime[ir:(ir+200)])) - meandt(precdt[ir:(ir+200)])
    buildup = range(start=0, length=length(precdt), step=DVM(delta)/(ir+600))
    vntime = @. precdt + Millisecond(round(Int64,buildup))
    return vntime # DateTime
end

"shift, interp, and stretch a time axis to match the GPS"
gpstime2gpsvndt(time, gpstime) = stretch_vn_to_gps( precise_dt(time, gpstime), gpstime )
"subtract leap seconds to get UTC time from GPS time"
gps2utc(gpsdt) = gpsdt - Second(leap_seconds(DateTime(gpsdt)))
gpstime2utcvndt(time, gpstime) = gps2utc( gpstime2gpsvndt(time, gpstime) )


## POSMV timing

"convert the GPS time from the POSMV PASHR message to a DateTime"
function pashr_gps_datetime(pashrtime, gpstime)
    gps_diff_12 = Millisecond(gpstime - (Time(pashrtime) - Hour(12)))
    gps_diff = mod(gps_diff_12, Day(1)) - Hour(12)
    gpsdt = pashrtime + gps_diff
    #DateTime(Date(pashrtime)) 
       #  + ( Millisecond.(gpstime - Time(pashrtime))
end



## functions to find the lag to maximize the autocorrelation between two time series

# e.g. the POSMV and VectorNav time series
# where VectorNav is factor=10 times the sampling rate of the POSMV. 
# Should also support factor=1.

# StatsBase.crosscor method for finding the optimal lag
# when time series are at same sampling rate.
"lag that maximizes the cross correlation"
function lag_xcor(x, y) 
    z = StatsBase.crosscor(x, y)
    delay=argmax(z)
    return lag = delay - div((size(z,1)+1),2) # lag=0 is in the center
end

"Cross-correlation function"
# use StateBase.crosscor
# could use DSP's trickier crosscor(x,y,lags), maybe faster using FFTs
function cross_correlation(lag, high_freq_series, low_freq_series, factor, start_idx, epoch_length)
    # start_idx is a low frequency index for the start of the epoch
    "index the epoch as a low frequency view"
    epx(lf) = lf[start_idx .+ (0:epoch_length-1)]
    # epx(lf) = lf[start_idx : min(start_idx+epoch_length-1, length(lf))] # ensure right in bounds

    # downsample and subset views
    if lag >= 0 # nonnegative lag
        # lag and downsample the high freq series
        x = epx( high_freq_series[round(Integer,lag+1):factor:end] )
        y = epx( low_freq_series )
    else # negative lag: indices need to be shifted in bounds
        # inbounds indices are >=1;
        # neg lag shifts both high and low freq series
        high_freq_offset = factor + lag%factor
        hfi = round(Integer, high_freq_offset%factor + 1)
        x = epx( high_freq_series[hfi:factor:end] ) # this line errors out of bounds!
        # positive offset for start of y
        low_freq_offset = round( Integer, -floor(lag / factor) )
        lfi = round(Integer, low_freq_offset + 1)
        y = epx( low_freq_series[lfi:end] ) # neg lag also shifts y
        # maxlen = max(length(y0), length(x0))
        # x = x0[1:maxlen] # shorten arrays if necessary
        # y = y0[1:maxlen]
    end
    cor(x, y) # StatsBase
end

"find the lag that maximizes the cross-correlation"
#=
function find_optimal_lag(f, a, b) # a and b are left and right bounds
    result = Optim.optimize(x -> -f(x), a, b, Optim.Brent()) # optim for max correlation
    # result = Optim.optimize(x -> -f(x), a, b, Optim.GoldenSection()) # either Brent or GoldenSection works
    return Optim.minimizer(result) # returns the argument of minimum of f
end
=#

function grid_search(high_freq_series, low_freq_series, factor, start_idx, epoch_length, initial_guess, window=40)   
    # find lag of max correlation with method from Optim, or brute force argmax
    crosscor_f(lag) = cross_correlation(lag, high_freq_series, low_freq_series, factor, start_idx, epoch_length)
    # find max correlation by brute force
    bruteforce_optimal_lag( lags::AbstractVector ) = lags[ argmax(crosscor_f.(lags)) ]
    "index the epoch from/as a low frequency view"
    # epx(lf) = lf[start_idx .+ (0:epoch_length-1)]
    eplen(lf) = min(epoch_length, length(lf)-start_idx+1) # emulates length of epx() windows
    
    # Initialize variables for the best lag and correlation
    left = initial_guess - window
    right = initial_guess + window
    # (inbounding was handed within cross_correlation)
    
    # check lagged (downsampled) views are in bounds
    nx = epoch_length
    ny = epoch_length
    for lag = [left, right]
        if lag >= 0 # nonnegative lag
            # lag and downsample the high freq series
            nx = min(nx, eplen( round(Integer,lag+1):factor:length(high_freq_series) ))
            ny = min(ny, eplen( 1:length(low_freq_series) ))
        else # negative lag: indices need to be shifted in bounds
            # inbounds indices are >=1;
            # neg lag shifts both high and low freq series
            high_freq_offset = factor + lag%factor
            hfi = round(Integer, high_freq_offset%factor + 1)
            # hfi(lag, factor=factor) = round(Integer, start_idx + (lag%factor + factor)%factor)  #nope not in one step
            nx = min(nx, eplen( hfi:factor:length(high_freq_series) )) # this line errors out of bounds!
            # positive offset for start of y
            low_freq_offset = round( Integer, -floor(lag / factor) )
            lfi = round(Integer, low_freq_offset + 1)
            ny = min(ny, eplen( lfi:length(low_freq_series) )) # neg lag also shifts y
        end
    end

    # condition on being inbounds
    if nx==ny==epoch_length # there is data to fill the epoch
        # optimal_lag = find_optimal_lag(crosscor_f, left, right) # doesn't work well
        optimal_lag = bruteforce_optimal_lag(left:right)
    elseif min(nx,ny) < epoch_length # epoch window out of range
        optimal_lag = missing
    end
    return optimal_lag #, best_cor
end

# not used; use grid_search directly
"""
Find the lag for local epochs
"""
function find_lags_iterative(high_freq_series, low_freq_series, factor, start_idxs, epoch_length, initial_guess; window=40)
    # optimal_lags = missing .+ zeros(Union{Missing, Int64}, size(start_idxs))
    optimal_lags = Vector{Union{Missing, Int64}}(missing, size(start_idxs))
    
    ic = 0
    for start_idx in start_idxs # loop each epoch
        if start_idx < length(low_freq_series) - epoch_length
            ic += 1
            optlag = grid_search(high_freq_series, low_freq_series, factor, start_idx, epoch_length, initial_guess, window)
            # Use the optimal lag from the current epoch as the 
            # initial guess for the next epoch

            if !ismissing(optlag)
                initial_guess = optlag
                optimal_lags[ic] = round(Int64, optlag)
            end
        end
    end
    
    return optimal_lags[1:ic]
end


end # module timing_lidar