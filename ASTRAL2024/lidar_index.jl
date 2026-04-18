module lidar_index

using Dates
using JLD2 

export LidarIndex, build_lidar_index, chunk_files, chunk_local_beams, dt_to_chunk, file_to_chunks, file_idx_to_chunks

struct LidarIndex
    # one entry per chunk (ic = 1:nchunks)
    ibeam_start  ::Vector{Int}      # global beam index of chunk start (= ist[ic])
    ibeam_end    ::Vector{Int}      # global beam index of chunk end   (= ien[ic])
    dt_start     ::Vector{DateTime} # dtime[ibeam_start[ic]]
    dt_end       ::Vector{DateTime} # dtime[ibeam_end[ic]]
    ifile_first  ::Vector{Int}      # index of first file spanning chunk
    ifile_last   ::Vector{Int}      # index of last file spanning chunk (= ifile_first if wholly within one file)

    # one entry per file (ifile = 1:nfiles)
    filenames          ::Vector{String}
    file_ibeam_start   ::Vector{Int}   # = bigind_file_start
    file_ibeam_end     ::Vector{Int}   # = bigind_file_end
    file_ichunk_first  ::Vector{Int}   # first chunk whose start falls in this file
    file_ichunk_last   ::Vector{Int}   # last chunk whose end falls in this file
end

# constructor
function build_lidar_index(dtime, ist, ien, files, file_ibeam_start, file_ibeam_end)
    nchunks = length(ist)
    nfiles  = length(files)

    ifile_first = [ findfirst(file_ibeam_start .<= ist[ic] .<= file_ibeam_end) for ic in 1:nchunks ]
    ifile_last  = [ findfirst(file_ibeam_start .<= ien[ic] .<= file_ibeam_end) for ic in 1:nchunks ]

    file_ichunk_first = [ findfirst(ifile_first .<= ifile) for ifile in 1:nfiles ]  # or 0 if none
    file_ichunk_last  = [ findlast( ifile_last  .>= ifile) for ifile in 1:nfiles ]

    LidarIndex(
        ist, ien,
        dtime[ist], dtime[ien],
        ifile_first, ifile_last,
        files,
        file_ibeam_start, file_ibeam_end,
        something.(file_ichunk_first, 0),
        something.(file_ichunk_last,  0)
    )
end

# lookup functions
# chunk ic → filename(s)
chunk_files(idx, ic) = idx.filenames[idx.ifile_first[ic]:idx.ifile_last[ic]]

# chunk ic → within-file beam range in its first file
chunk_local_beams(idx, ic) = let f = idx.ifile_first[ic]
    (idx.ibeam_start[ic] - idx.file_ibeam_start[f] + 1,
     idx.ibeam_end[ic]   - idx.file_ibeam_start[f] + 1)   # clipped at file boundary if spans
end

# datetime → chunk index
dt_to_chunk(idx, dt) = findlast(idx.dt_start .<= dt)

# filename → chunk indices
file_to_chunks(idx, fname) = let ifile = findfirst(==(fname), idx.filenames)
    idx.file_ichunk_first[ifile]:idx.file_ichunk_last[ifile]
end

# file index → chunk indices
file_idx_to_chunks(idx, ifile) = idx.file_ichunk_first[ifile]:idx.file_ichunk_last[ifile]

end # module lidar_index