#============================================================#
##################### Reading PottsGraph #####################
#============================================================#

# PottsGraph(file::AbstractString) = read_graph_extended(file)
"""
    read_graph
    read_potts_graph

Read `PottsGraph` object from a file.
"""
read_graph(file::AbstractString) = read_graph_extended(file)
"""
    read_potts_graph

Alias for `read_graph`.
"""
read_potts_graph = read_graph

function read_graph_extended(file)
    ## Go through file twice: first to get L and q, the second to store parameters
    q = 0
    L = 0
    min_idx = Inf
    index_style = 1
    for line in eachline(file)
        @assert is_valid_line(line) "Format problem with line:\n $line \n--> Expected `J i j a b` or `h i a`."
        if !isempty(line) && line[1] == 'h'
            i, a, val = parse_field_line(line)
            if i > L
                L = i
            end
            if a > q
                q = a
            end
            if i < min_idx || a < min_idx
                min_idx = min(i, a)
            end
        end
    end
    @assert min_idx == 1 || min_idx == 0 "Issue with indexing: smallest index found is $min_idx"
    index_style = (min_idx == 0 ? 0 : 1)
    if index_style == 0
        L += 1
        q += 1
    end

    g = PottsGraph(L, q)
    for line in eachline(file)
        if line[1] == 'J'
            i, j, a, b, val = parse_coupling_line(line)
            index_style == 0 && (i += 1; j += 1; a += 1; b += 1)
            g.J[a,b,i,j] = val
            g.J[b,a,j,i] = val
        elseif line[1] == 'h'
            i, a, val = parse_field_line(line)
            index_style == 0 && (i += 1; a += 1)
            g.h[a,i] = val
        end
    end

    return g
end

function is_valid_line(line)
    if isnothing(match(r"J [0-9]+ [0-9]+ [0-9]+ [0-9]+", line)) &&
        isnothing(match(r"h [0-9]+ [0-9]+", line))
        return false
    else
        return true
    end
end
function parse_field_line(line)
    s = split(line, " ")
    i = parse(Int, s[2])
    a = parse(Int, s[3])
    val = parse(Float64, s[4])
    return i, a, val
end
function parse_coupling_line(line)
    s = split(line, " ")
    i = parse(Int, s[2])
    j = parse(Int, s[3])
    a = parse(Int, s[4])
    b = parse(Int, s[5])
    val = parse(Float64, s[6])
    return i, j, a, b, val
end

#============================================================#
##################### Writing PottsGraph #####################
#============================================================#

"""
    write(file::AbstractString, g::PottsGraph; sigdigits)

Write parameters of `g` to `file` using the format `J i j a b value`.
"""
function write(file::AbstractString, g::PottsGraph; sigdigits=5, index_style=1)
    return write_graph_extended(file, g, sigdigits, index_style)
end

function write_graph_extended(file::AbstractString, g::PottsGraph, sigdigits, index_style)
    @assert index_style == 0 || index_style == 1 "Got `index_style==`$(index_style)"
    L, q = size(g)
    open(file, "w") do f
        for i in 1:L, j in (i+1):L, a in 1:q, b in 1:q
            val = round(g.J[a,b,i,j]; sigdigits)
            if index_style == 0
                write(f, "J $(i-1) $(j-1) $(a-1) $(b-1) $val\n")
            elseif index_style == 1
                write(f, "J $i $j $a $b $val\n")
            end
        end
        for i in 1:L, a in 1:q
            val = round(g.h[a,i]; sigdigits)
            if index_style == 0
                write(f, "h $(i-1) $(a-1) $val\n" )
            elseif index_style == 1
                write(f, "h $i $a $val\n" )
            end
        end
    end
end
